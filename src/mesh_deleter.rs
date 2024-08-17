use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    ops::{Index, IndexMut},
    process::id,
};

use linear_isomorphic::RealField;
use num::Float;

use crate::{
    container_trait::{PrimitiveContainer, RedgeContainers, VertData},
    face_handle,
    helpers::{
        disable_edge_meta, disable_face_meta, disable_hedge_meta, disable_vert_meta,
        join_radial_cycles, join_vertex_cycles, pick_different_edge, remove_edge_from_cycle,
        remove_hedge_from_face, remove_hedge_from_radial,
    },
    validation::{correctness_state, RedgeCorrectness},
    wavefront_loader::ObjData,
    EdgeId, FaceId, HedgeId, Redge, StarCycleNode, VertId,
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
    mesh: Redge<R>,
}

impl<R: RedgeContainers> MeshDeleter<R> {
    pub fn start_deletion(mesh: Redge<R>) -> Self {
        Self { mesh }
    }

    pub fn end_deletion(mut self) -> Redge<R> {
        let (vert_frag, edge_frag, hedge_frag, face_frag) = self.compute_fragmentation_maps();
        self.apply_defragmentation_maps(vert_frag, edge_frag, hedge_frag, face_frag);
        self.mesh
    }

    pub fn mesh(&mut self) -> &mut Redge<R> {
        &mut self.mesh
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

        let mut edge_fragmentation = HashMap::<EdgeId, usize>::new();
        let mut counter = 0;
        for e in self.mesh.edges_meta.iter() {
            if e.is_active {
                edge_fragmentation.insert(e.id, counter);
                counter += 1;
            }
        }

        let mut hedge_fragmentation = HashMap::<HedgeId, usize>::new();
        let mut counter = 0;
        for h in self.mesh.hedges_meta.iter() {
            if h.is_active {
                hedge_fragmentation.insert(h.id, counter);
                counter += 1;
            }
        }

        let mut face_fragmentation = HashMap::<FaceId, usize>::new();
        let mut counter = 0;
        for f in self.mesh.faces_meta.iter() {
            if f.is_active {
                face_fragmentation.insert(f.id, counter);
                counter += 1;
            }
        }

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

            if !e.hedge_id.is_absent() {
                e.hedge_id = HedgeId(*hedge_frag.get(&e.hedge_id).unwrap());
            }

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
    }

    pub fn remove_edge(&mut self, edge_id: EdgeId) {
        let incident_faces: Vec<_> = self
            .mesh
            .edge_handle(edge_id)
            .hedge()
            .radial_neighbours()
            .map(|h| h.face().id())
            .collect();

        let handle = self.mesh.edge_handle(edge_id);

        let [v1, v2] = handle.vertex_ids();
        let v1_safe_edge = handle
            .v1()
            .star_edges()
            .find(|e| e.id() != edge_id)
            .unwrap()
            .id();

        let v2_safe_edge = handle
            .v2()
            .star_edges()
            .find(|e| e.id() != edge_id)
            .unwrap()
            .id();

        for fid in incident_faces {
            self.remove_face(fid);
        }

        remove_edge_from_cycle(edge_id, crate::Endpoint::V1, &mut self.mesh);
        remove_edge_from_cycle(edge_id, crate::Endpoint::V2, &mut self.mesh);

        disable_edge_meta(edge_id, &mut self.mesh);

        self.mesh.verts_meta[v1.to_index()].edge_id = v1_safe_edge;
        self.mesh.verts_meta[v2.to_index()].edge_id = v2_safe_edge;
    }

    // Note: There's no doubt this could be made more efficient, but
    // protecting the invariants is very hard. Don't touch this function
    // unless there's a REALLY compelling case it needs to be done.
    pub fn collapse_edge(&mut self, edge_id: EdgeId) -> VertId {
        let v1 = self.mesh.edge_handle(edge_id).v1().id();
        let v2 = self.mesh.edge_handle(edge_id).v2().id();

        let edges_to_merge: Vec<_> = self
            .mesh
            .edge_handle(edge_id)
            .hedge()
            .radial_neighbours()
            .map(|h| [h.face_prev().edge().id(), h.face_next().edge().id()])
            .collect();

        // TODO: This is incorrect for non-triangular faces as one does not need to
        // remove the faces. Additional logic is needed to keep or restore those faces.
        self.remove_edge(edge_id);
        debug_assert!(correctness_state(&self.mesh) == RedgeCorrectness::Correct,);

        let v2_edges: Vec<_> = self
            .mesh
            .vert_handle(v2)
            .star_edges()
            .map(|e| e.id())
            .collect();

        let v1_edges: Vec<_> = self
            .mesh
            .vert_handle(v1)
            .star_edges()
            .map(|e| e.id())
            .collect();

        // Verify that connectivity is sane.
        debug_assert!({
            let mut status = true;
            for e in &v2_edges {
                let edge = &self.mesh.edges_meta[e.to_index()];
                if !edge.is_active {
                    continue;
                }
                status = status && (edge.vert_ids[0] == v2 || edge.vert_ids[1] == v2);
            }
            status
        });

        // Make all things touching V2 now touch v1.
        for eid in &v2_edges {
            let radials: Vec<_> = self
                .mesh
                .edge_handle(*eid)
                .hedge()
                .radial_neighbours()
                .map(|h| h.id())
                .collect();
            for radial in radials {
                if self.mesh.hedges_meta[radial.to_index()].source_id == v2 {
                    self.mesh.hedges_meta[radial.to_index()].source_id = v1;
                }
            }
            *self.mesh.edges_meta[eid.to_index()].at(v2) = v1;
        }

        // Merge the endpoint vertex cycle.
        let edges: Vec<_> = v2_edges.iter().chain(v1_edges.iter()).collect();

        for i in 0..=edges.len() {
            let ep = edges[i % edges.len()];
            let en = edges[(i + 1) % edges.len()];

            self.mesh.edges_meta[ep.to_index()].cycle_mut(v1).next_edge = *en;
            self.mesh.edges_meta[en.to_index()].cycle_mut(v1).prev_edge = *ep;

            self.mesh.verts_meta[v1.to_index()].edge_id = *ep;
        }

        disable_vert_meta(v2, &mut self.mesh);

        for [e1, e2] in edges_to_merge {
            let h1 = self.mesh.edges_meta[e1.to_index()].hedge_id;
            let h2 = self.mesh.edges_meta[e2.to_index()].hedge_id;

            let [e2_v1, e2_v2] = self.mesh.edges_meta[e1.to_index()].vert_ids;

            // Make all hedges of e2 point to e1.
            let orbit: Vec<_> = self
                .mesh
                .edge_handle(e2)
                .hedge()
                .radial_neighbours()
                .map(|h| h.id())
                .collect();
            for h in orbit {
                self.mesh.hedges_meta[h.to_index()].edge_id = e1;
            }

            join_radial_cycles(h1, h2, &mut self.mesh);
            remove_edge_from_cycle(e2, crate::Endpoint::V1, &mut self.mesh);
            remove_edge_from_cycle(e2, crate::Endpoint::V2, &mut self.mesh);
            disable_edge_meta(e2, &mut self.mesh);

            // Prevent invalidation of the endpoints if they pointed to the removed edge.
            self.mesh.verts_meta[e2_v1.to_index()].edge_id = e1;
            self.mesh.verts_meta[e2_v2.to_index()].edge_id = e1;
        }

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
        // deleter.remove_edge(EdgeId(0));
        deleter.collapse_edge(EdgeId(0));

        let state = manifold_state(deleter.mesh());
        debug_assert!(state == RedgeManifoldness::IsManifold, "{:?}", state);

        let redge = deleter.end_deletion();

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/loop_cube.obj");
    }
}
