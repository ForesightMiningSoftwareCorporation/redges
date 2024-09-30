use std::{
    collections::HashMap,
    ops::{Index, Mul},
    os::linux::raw::stat,
    usize,
};

use linear_isomorphic::{ArithmeticType, InnerSpace, RealField};
use num_traits::float::TotalOrder;
use ordered_float::FloatCore;

use crate::{
    container_trait::{PrimitiveContainer, RedgeContainers, VertData},
    edge_handle,
    face_handle::FaceMetrics,
    helpers::{
        disable_edge_meta, disable_face_meta, disable_hedge_meta, disable_vert_meta,
        fix_digon_face, hedge_collapse, join_radial_cycles, join_vertex_cycles,
        join_vertex_cycles_at, remove_edge_from_cycle, remove_hedge_from_radial,
    },
    validation::{correctness_state, RedgeCorrectness},
    wavefront_loader::ObjData,
    EdgeId, Endpoint, FaceId, HedgeId, Redge, VertId,
};

/// Helper macro to avoid redundancy when defragmenting mesh data.
macro_rules! compact_mesh_data {
    ($self:ident, $data:ident, $meta:ident) => {
        debug_assert!(
            $self.mesh.$data.len() == $self.mesh.$meta.len() || $self.mesh.$data.len() == 0
        );
        let mut count = 0;
        for i in 0..$self.mesh.$meta.len() {
            if !$self.mesh.$meta[i].is_active {
                continue;
            };

            $self.mesh.$meta[count] = $self.mesh.$meta[i].clone();
            let data = $self.mesh.$data.get(i as u64).clone();

            if $self.mesh.$data.len() > 0 {
                *$self.mesh.$data.get_mut(count as u64) = data;
            }

            count += 1;
        }

        $self.mesh.$meta = $self.mesh.$meta[0..count].to_vec();
        $self.mesh.$data.resize(count.min($self.mesh.$data.len()));
    };
}

pub struct MeshDeleter<R: RedgeContainers> {
    pub(crate) mesh: Redge<R>,
    deleted_verts: usize,
    deleted_edges: usize,
    deleted_faces: usize,
}

impl<R: RedgeContainers> MeshDeleter<R> {
    pub fn start_deletion(mesh: Redge<R>) -> Self {
        Self {
            mesh,
            deleted_verts: 0,
            deleted_edges: 0,
            deleted_faces: 0,
        }
    }

    pub fn end_deletion(mut self) -> Redge<R> {
        let (vert_frag, edge_frag, hedge_frag, face_frag) = self.compute_fragmentation_maps();
        self.apply_defragmentation_maps(vert_frag, edge_frag, hedge_frag, face_frag);
        self.mesh
    }

    pub fn mesh(&mut self) -> &mut Redge<R> {
        &mut self.mesh
    }

    pub fn active_vert_count(&self) -> usize {
        self.mesh.vert_count() - self.deleted_verts
    }

    pub fn active_edge_count(&self) -> usize {
        self.mesh.edge_count() - self.deleted_edges
    }

    pub fn active_face_count(&self) -> usize {
        self.mesh.face_count() - self.deleted_faces
    }

    pub fn compute_fragmentation_maps(
        &self,
    ) -> (
        HashMap<VertId, usize>,
        HashMap<EdgeId, usize>,
        HashMap<HedgeId, usize>,
        HashMap<FaceId, usize>,
    ) {
        let mut vertex_fragmentation = HashMap::<VertId, usize>::new();
        let mut counter = 0;
        for v in self.mesh.verts_meta.iter() {
            if v.is_active {
                vertex_fragmentation.insert(v.id, counter);
                counter += 1;
            }
        }
        vertex_fragmentation.insert(VertId::ABSENT, usize::MAX);

        let mut edge_fragmentation = HashMap::<EdgeId, usize>::new();
        let mut counter = 0;
        for e in self.mesh.edges_meta.iter() {
            if e.is_active {
                edge_fragmentation.insert(e.id, counter);
                counter += 1;
            }
        }
        edge_fragmentation.insert(EdgeId::ABSENT, usize::MAX);

        let mut hedge_fragmentation = HashMap::<HedgeId, usize>::new();
        let mut counter = 0;
        for h in self.mesh.hedges_meta.iter() {
            if h.is_active {
                hedge_fragmentation.insert(h.id, counter);
                counter += 1;
            }
        }
        hedge_fragmentation.insert(HedgeId::ABSENT, usize::MAX);

        let mut face_fragmentation = HashMap::<FaceId, usize>::new();
        let mut counter = 0;
        for f in self.mesh.faces_meta.iter() {
            if f.is_active {
                face_fragmentation.insert(f.id, counter);
                counter += 1;
            }
        }
        face_fragmentation.insert(FaceId::ABSENT, usize::MAX);

        (
            vertex_fragmentation,
            edge_fragmentation,
            hedge_fragmentation,
            face_fragmentation,
        )
    }

    fn apply_defragmentation_maps(
        &mut self,
        vert_frag: HashMap<VertId, usize>,
        edge_frag: HashMap<EdgeId, usize>,
        hedge_frag: HashMap<HedgeId, usize>,
        face_frag: HashMap<FaceId, usize>,
    ) {
        for v in self.mesh.verts_meta.iter_mut() {
            if !v.is_active {
                continue;
            }
            v.id = VertId(*vert_frag.get(&v.id).unwrap());
            v.edge_id = EdgeId(*edge_frag.get(&v.edge_id).unwrap());
        }

        for e in self.mesh.edges_meta.iter_mut() {
            if !e.is_active {
                continue;
            }
            e.id = EdgeId(*edge_frag.get(&e.id).unwrap());

            e.vert_ids[0] = VertId(*vert_frag.get(&e.vert_ids[0]).unwrap());
            e.vert_ids[1] = VertId(*vert_frag.get(&e.vert_ids[1]).unwrap());

            e.hedge_id = HedgeId(*hedge_frag.get(&e.hedge_id).unwrap());

            e.v1_cycle.next_edge = EdgeId(*edge_frag.get(&e.v1_cycle.next_edge).unwrap());
            e.v1_cycle.prev_edge = EdgeId(*edge_frag.get(&e.v1_cycle.prev_edge).unwrap());

            e.v2_cycle.next_edge = EdgeId(*edge_frag.get(&e.v2_cycle.next_edge).unwrap());
            e.v2_cycle.prev_edge = EdgeId(*edge_frag.get(&e.v2_cycle.prev_edge).unwrap());
        }

        for h in self.mesh.hedges_meta.iter_mut() {
            if !h.is_active {
                continue;
            }

            h.id = HedgeId(*hedge_frag.get(&h.id).unwrap());
            h.source_id = VertId(*vert_frag.get(&h.source_id).unwrap());
            h.edge_id = EdgeId(*edge_frag.get(&h.edge_id).unwrap());
            h.radial_next_id = HedgeId(*hedge_frag.get(&h.radial_next_id).unwrap());
            h.radial_prev_id = HedgeId(*hedge_frag.get(&h.radial_prev_id).unwrap());
            h.face_next_id = HedgeId(*hedge_frag.get(&h.face_next_id).unwrap());
            h.face_prev_id = HedgeId(*hedge_frag.get(&h.face_prev_id).unwrap());
            h.face_id = FaceId(*face_frag.get(&h.face_id).unwrap());
        }

        for f in self.mesh.faces_meta.iter_mut() {
            if !f.is_active {
                continue;
            }

            f.id = FaceId(*face_frag.get(&f.id).unwrap());
            f.hedge_id = HedgeId(*hedge_frag.get(&f.hedge_id).unwrap());
        }

        compact_mesh_data!(self, vert_data, verts_meta);
        compact_mesh_data!(self, edge_data, edges_meta);
        compact_mesh_data!(self, face_data, faces_meta);

        let mut count = 0;
        for i in 0..self.mesh.hedges_meta.len() {
            if !self.mesh.hedges_meta[i].is_active {
                continue;
            };

            self.mesh.hedges_meta[count] = self.mesh.hedges_meta[i].clone();

            count += 1;
        }

        self.mesh.hedges_meta.truncate(count);
    }

    pub fn remove_face(&mut self, face_id: FaceId) {
        let face_hedges: Vec<_> = self
            .mesh
            .face_handle(face_id)
            .hedge()
            .face_loop()
            .map(|h| h.id())
            .collect();

        for hid in &face_hedges {
            let p = self.mesh.hedges_meta[hid.to_index()].radial_next_id;
            let e = self.mesh.hedges_meta[p.to_index()].edge_id;
            remove_hedge_from_radial(*hid, &mut self.mesh);

            debug_assert!(self.mesh.hedges_meta[p.to_index()].edge_id == e);
        }

        for hid in face_hedges {
            disable_hedge_meta(hid, &mut self.mesh);
        }

        disable_face_meta(face_id, &mut self.mesh);

        self.deleted_faces += 1;
    }

    pub fn remove_edge(&mut self, edge_id: EdgeId) {
        let edge_handle = self.mesh.edge_handle(edge_id);
        let incident_faces: Vec<_> = if edge_handle.has_hedge() {
            edge_handle
                .hedge()
                .radial_neighbours()
                .map(|h| h.face().id())
                .collect()
        } else {
            Vec::new()
        };

        let [v1, v2] = edge_handle.vertex_ids();
        let v1_safe_edge = edge_handle
            .v1()
            .star_edges()
            .find(|e| e.id() != edge_id)
            .map_or(EdgeId::ABSENT, |e| e.id());

        let v2_safe_edge = edge_handle
            .v2()
            .star_edges()
            .find(|e| e.id() != edge_id)
            .map_or(EdgeId::ABSENT, |e| e.id());

        for fid in incident_faces {
            self.remove_face(fid);
        }

        remove_edge_from_cycle(edge_id, crate::Endpoint::V1, &mut self.mesh);
        remove_edge_from_cycle(edge_id, crate::Endpoint::V2, &mut self.mesh);

        disable_edge_meta(edge_id, &mut self.mesh);

        self.mesh.verts_meta[v1.to_index()].edge_id = v1_safe_edge;
        self.mesh.verts_meta[v2.to_index()].edge_id = v2_safe_edge;

        if v1_safe_edge == EdgeId::ABSENT {
            disable_vert_meta(v1, &mut self.mesh);
        }

        if v2_safe_edge == EdgeId::ABSENT {
            disable_vert_meta(v2, &mut self.mesh);
        }

        self.deleted_edges += 1;
    }

    // Note: There's no doubt this could be made more efficient, but
    // protecting the invariants is very hard. Don't touch this function
    // unless there's a REALLY compelling case it needs to be done.
    //
    // If you plan on modifying this function please read `docs/redge.pdf`.
    //
    /// Currently only works for triangular faces. Only call on edges that
    /// have faces pointing to them.
    pub fn collapse_edge<S>(&mut self, edge_id: EdgeId) -> VertId
    where
        VertData<R>: Index<usize, Output = S>,
        S: RealField,
    {
        // Collect all necessary elements before breaking the topology.
        let edge_handle = self.mesh.edge_handle(edge_id);
        let hedges_to_collapse: Vec<_> = edge_handle
            .hedge()
            .radial_neighbours()
            .map(|f| f.id())
            .collect();

        let v1_edges: Vec<_> = edge_handle
            .v1()
            .star_edges()
            .map(|e| e.id())
            .filter(|id| *id != edge_id)
            .collect();
        let v2_edges: Vec<_> = edge_handle
            .v2()
            .star_edges()
            .map(|e| e.id())
            .filter(|id| *id != edge_id)
            .collect();

        let v1 = edge_handle.v1().id();
        let v2 = edge_handle.v2().id();

        // For safety, make sure that v1 does not point to the current edge. Since we are about to remove it.
        self.mesh.verts_meta[v1.to_index()].edge_id = v1_edges[0];

        remove_edge_from_cycle(edge_id, Endpoint::V1, &mut self.mesh);
        remove_edge_from_cycle(edge_id, Endpoint::V2, &mut self.mesh);

        let mut faces = Vec::new();
        // Join the hedges of each face.
        for hid in hedges_to_collapse {
            faces.push(self.mesh.hedges_meta[hid.to_index()].face_id);
            hedge_collapse(hid, &mut self.mesh);
        }

        // Update the edges incident on v2 to point to v1 instead.
        for eid in &v2_edges {
            *self.mesh.edges_meta[eid.to_index()].at(v2) = v1;

            let hedges: Vec<HedgeId> = self
                .mesh
                .edge_handle(*eid)
                .hedge()
                .radial_neighbours()
                .map(|h| h.id())
                .collect();

            // Make sure that any hedges orbitting this edge have their sources updated appropriately.
            for hid in hedges {
                if self.mesh.hedges_meta[hid.to_index()].source_id == v2 {
                    self.mesh.hedges_meta[hid.to_index()].source_id = v1;
                }
            }
        }

        join_vertex_cycles_at(v1_edges[0], v2_edges[0], v1, &mut self.mesh);

        disable_edge_meta(edge_id, &mut self.mesh);
        disable_vert_meta(v2, &mut self.mesh);

        for face in faces {
            if self.mesh.face_handle(face).side_count() == 2 {
                fix_digon_face(face, &mut self.mesh);
                self.deleted_faces += 1;
                self.deleted_edges += 1;
            }
        }

        debug_assert!(correctness_state(&self.mesh) == RedgeCorrectness::Correct);

        v1
    }
}

#[cfg(test)]
mod tests {
    use crate::validation::{manifold_state, RedgeManifoldness};
    use crate::wavefront_loader::ObjData;

    use super::*;

    #[test]
    fn test_edge_collapse() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/loop_cube.obj");

        let redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let state = manifold_state(&redge);
        debug_assert!(state == RedgeManifoldness::IsManifold, "{:?}", state);

        let mut deleter = MeshDeleter::start_deletion(redge);
        deleter.collapse_edge(EdgeId(0));

        let state = manifold_state(deleter.mesh());
        debug_assert!(state == RedgeManifoldness::IsManifold, "{:?}", state);

        let redge = deleter.end_deletion();

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/loop_cube.obj");
    }

    #[test]
    fn test_edge_collapse_experimental() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/loop_cube.obj");

        let redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let state = manifold_state(&redge);
        debug_assert!(state == RedgeManifoldness::IsManifold, "{:?}", state);

        let mut deleter = MeshDeleter::start_deletion(redge);
        deleter.collapse_edge(EdgeId(0));

        let state = manifold_state(deleter.mesh());
        debug_assert!(state == RedgeManifoldness::IsManifold, "{:?}", state);

        let redge = deleter.end_deletion();

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/loop_cube.obj");
    }
}
