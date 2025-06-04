//! Wrapper to enable destructive operations on a radial edge.
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    ops::Index,
};

use linear_isomorphic::RealField;

use crate::{
    container_trait::{
        FaceAttributeGetter, FaceData, PrimitiveContainer, RedgeContainers, VertData,
    },
    face_handle::FaceDegeneracies,
    helpers::{
        digon_holes_to_edge, digon_to_edge, disable_edge_meta, disable_face_meta,
        disable_hedge_meta, disable_vert_meta, hedge_collapse, join_vertex_cycles_at,
        remove_edge_from_cycle, remove_hedge_from_radial,
    },
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

/// Wrapper around a radial edge, enabling destructive operations.
pub struct MeshDeleter<R: RedgeContainers> {
    pub(crate) mesh: Redge<R>,
    deleted_verts: usize,
    deleted_edges: usize,
    deleted_faces: usize,
}

type FragmentationMaps = (
    HashMap<VertId, usize>,
    HashMap<EdgeId, usize>,
    HashMap<HedgeId, usize>,
    HashMap<FaceId, usize>,
);
impl<R: RedgeContainers> MeshDeleter<R> {
    /// Construct a deleter and hold the radial edge hostage until done.
    pub fn start_deletion(mesh: Redge<R>) -> Self {
        Self {
            mesh,
            deleted_verts: 0,
            deleted_edges: 0,
            deleted_faces: 0,
        }
    }

    /// Finish deleting and release the hostage radial edge.
    pub fn end_deletion(mut self) -> Redge<R> {
        let (vert_frag, edge_frag, hedge_frag, face_frag) = self.compute_fragmentation_maps();
        self.apply_defragmentation_maps(vert_frag, edge_frag, hedge_frag, face_frag);
        self.mesh
    }
    /// Get the underlying radial edge, use with care.
    pub fn mesh(&mut self) -> &mut Redge<R> {
        &mut self.mesh
    }
    /// Count of *actual* active vertices (i.e. after removing deleted ones).
    pub fn active_vert_count(&self) -> usize {
        self.mesh.vert_count() - self.deleted_verts
    }
    /// Count of *actual* active edges (i.e. after removing deleted ones).
    pub fn active_edge_count(&self) -> usize {
        self.mesh.edge_count() - self.deleted_edges
    }
    /// Count of *actual* active faces (i.e. after removing deleted ones).
    pub fn active_face_count(&self) -> usize {
        self.mesh.face_count() - self.deleted_faces
    }
    /// Maps between old ids and ids after applying defragmentation.
    pub fn compute_fragmentation_maps(&self) -> FragmentationMaps {
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

    /// Remove a face from the underlying radial edge.
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

    /// Remove an edge from the underlying radial edge.
    pub fn remove_edge(&mut self, edge_id: EdgeId) {
        let edge_handle = self.mesh.edge_handle(edge_id);
        let incident_faces: Vec<_> = if edge_handle.has_hedge() {
            edge_handle
                .hedge()
                .radial_loop()
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

    /// Remove a vertex from the underlying radial edge.
    pub fn remove_vert(&mut self, vert_id: VertId) {
        let faces_to_remove: Vec<_> = self
            .mesh
            .vert_handle(vert_id)
            .incident_faces()
            .map(|f| f.id())
            .collect();

        for fid in faces_to_remove {
            self.remove_face(fid);
        }

        let edges_to_remove: Vec<_> = self
            .mesh
            .vert_handle(vert_id)
            .star_edges()
            .map(|f| f.id())
            .collect();

        for eid in edges_to_remove {
            self.remove_edge(eid);
        }

        self.deleted_verts += 1;
        disable_vert_meta(vert_id, &mut self.mesh);
    }

    /// Inspect the mesh for overlapping faces and remove them. This is
    /// expensive so only use it if you can't demonstrate that the mesh is manifold.
    pub fn remove_overlapping_faces(&mut self) {
        let faces_ids = self.mesh.meta_faces().map(|f| f.id()).collect::<Vec<_>>();

        for fid in faces_ids {
            let handle = self.mesh.face_handle(fid);
            if !handle.is_active() {
                continue;
            }

            let status = handle.check_degeneracies();
            if status == FaceDegeneracies::Doppelganger {
                self.remove_face(fid);
            }
        }
    }

    // Note: There's no doubt this could be made more efficient, but
    // protecting the invariants is very hard. Don't touch this function
    // unless there's a REALLY compelling case it needs to be done.
    //
    // If you plan on modifying this function please read `docs/redge.pdf`.
    //
    /// Collapse an edge.
    /// Warning: Currently only works for triangular faces.
    fn collapse_edge<S>(&mut self, edge_id: EdgeId) -> VertId
    where
        VertData<R>: Index<usize, Output = S>,
        S: RealField,
    {
        // Collect all necessary elements before breaking the topology.
        let edge_handle = self.mesh.edge_handle(edge_id);
        let hedges_to_collapse: Vec<_> =
            edge_handle.hedge().radial_loop().map(|h| h.id()).collect();

        let mut v1_edges: Vec<_> = edge_handle
            .v1()
            .star_edges()
            .map(|e| e.id())
            .filter(|id| *id != edge_id)
            .collect();
        let mut v2_edges: Vec<_> = edge_handle
            .v2()
            .star_edges()
            .map(|e| e.id())
            .filter(|id| *id != edge_id)
            .collect();

        // Find and delete any edges that also connect these two endpoints.
        let degenerate = v1_edges
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .intersection(&v2_edges.iter().copied().collect::<BTreeSet<_>>())
            .copied()
            .collect::<BTreeSet<_>>();

        for eid in &degenerate {
            self.remove_edge(*eid);
        }
        // Keep only sane edges.
        v1_edges.retain(|eid| !degenerate.contains(eid));
        v2_edges.retain(|eid| !degenerate.contains(eid));

        let edge_handle = self.mesh.edge_handle(edge_id);

        let v1 = edge_handle.v1().id();
        let v2 = edge_handle.v2().id();

        // For safety, make sure that v1 does not point to the current edge. Since we are about to remove it.
        self.mesh.verts_meta[v1.to_index()].edge_id = if !v1_edges.is_empty() {
            v1_edges[0]
        } else if !v2_edges.is_empty() {
            v2_edges[0]
        } else {
            EdgeId::ABSENT
        };

        remove_edge_from_cycle(edge_id, Endpoint::V1, &mut self.mesh);
        remove_edge_from_cycle(edge_id, Endpoint::V2, &mut self.mesh);

        // Join the hedges of each face.
        for hid in hedges_to_collapse {
            hedge_collapse(hid, &mut self.mesh);
        }

        // Update the edges incident on v2 to point to v1 instead.
        for eid in &v2_edges {
            *self.mesh.edges_meta[eid.to_index()].at(v2) = v1;

            let hedges: Vec<HedgeId> = self
                .mesh
                .edge_handle(*eid)
                .hedge()
                .radial_loop()
                .map(|h| h.id())
                .collect();

            // Make sure that any hedges orbitting this edge have their sources updated appropriately.
            for hid in hedges {
                if self.mesh.hedges_meta[hid.to_index()].source_id == v2 {
                    self.mesh.hedges_meta[hid.to_index()].source_id = v1;
                }
            }
        }

        if !v1_edges.is_empty() && !v2_edges.is_empty() {
            join_vertex_cycles_at(v1_edges[0], v2_edges[0], v1, &mut self.mesh);
        }

        disable_edge_meta(edge_id, &mut self.mesh);
        disable_vert_meta(v2, &mut self.mesh);
        self.deleted_verts += 1;

        v1
    }

    /// Collapse an edge, and fix any created degenerate faces.
    /// Warning: Currently only works for triangular faces.
    pub fn collapse_edge_and_fix<S>(&mut self, edge_id: EdgeId) -> VertId
    where
        VertData<R>: Index<usize, Output = S>,
        S: RealField,
        FaceData<R>: FaceAttributeGetter<S>,
    {
        // Collect all necessary elements before breaking the topology.
        let edge_handle = self.mesh.edge_handle(edge_id);
        let faces = edge_handle
            .hedge()
            .radial_loop()
            .map(|h| h.face().id())
            .collect::<Vec<_>>();

        let vid = self.collapse_edge(edge_id);

        // TODO: This should be the true logic, however, on the stanford dragon,
        // there is one particular vertex where edges keep getting collapsed
        // in a weird way that creates a vertex with a large valence.
        // This doesn't happen when we follow the original radial order of the faces.
        // This is likely indicative that something in the quadrics simplification
        // is very sensitive to order, but I am not sure what and time is a resource.
        // Keep this in mind, we will likely need to revisit this bug.

        // // Collapsing an edge touching triangular faces will create digons.
        // // find all of them and turn them into edges, for sanity.
        // let digons: Vec<_> = self
        //     .mesh
        //     .vert_handle(vid)
        //     .incident_faces()
        //     .filter(|f| f.side_count() == 2)
        //     .map(|f| f.id())
        //     .collect();

        for face in faces {
            if self.mesh.face_handle(face).side_count() == 2 {
                digon_to_edge(face, &mut self.mesh);
                self.deleted_faces += 1;
                self.deleted_edges += 1;
            }
        }

        // Collapsing an edge in a triangular hole creates a hole with only two edges.
        // These holes are invisible and really hard to deal with, so just close them.
        let mut digon_holes = BTreeMap::new();
        for edge in self
            .mesh
            .vert_handle(vid)
            .star_edges()
            .filter(|e| e.is_boundary())
        {
            // A digon hole is defined by two boundary edges with the same endpoints.
            let val = digon_holes.entry(edge.opposite(vid)).or_insert(vec![]);
            val.push(edge.id());
        }

        digon_holes.retain(|_, l| l.len() >= 2);

        for (_, digons) in digon_holes {
            digon_holes_to_edge(digons, &mut self.mesh);
        }

        self.clean_doppleganger_faces(vid);

        vid
    }

    fn clean_doppleganger_faces(&mut self, vid: VertId) {
        // Find all faces which overlap in space (i.e. identical faces).
        let dihedron_faces: Vec<_> = self
            .mesh
            .vert_handle(vid)
            .incident_faces()
            .filter(|f| f.check_degeneracies() == FaceDegeneracies::Doppelganger)
            .collect();

        let mut victim_dihedron_faces = BTreeSet::new();
        for dihedron in dihedron_faces {
            // See if we have already victimized the other face in the dihedron. If we have not
            // then we victimize this face.
            if !victim_dihedron_faces.contains(&dihedron.hedge().radial_next().face().id()) {
                victim_dihedron_faces.insert(dihedron.id());
            }
        }

        // Kill one of the dihedrons.
        for victim_dihedron in victim_dihedron_faces {
            self.remove_face(victim_dihedron);
        }
    }

    /// Re-associate faces with their new vertex ids.
    pub fn update_face_corners<S: RealField>(&mut self, vid: VertId)
    where
        FaceData<R>: FaceAttributeGetter<S>,
    {
        // Find all faces incident on this vertex.
        let fids: Vec<_> = self
            .mesh
            .vert_handle(vid)
            .incident_faces()
            .map(|f| f.id())
            .collect();

        // For each such face, see if we need to update one of its corners.
        for fid in fids {
            let face_handle = self.mesh.face_handle(fid);
            // Get the current topology ids, as well as the stored ids inside the face data.
            let vert_ids: Vec<_> = face_handle.vertex_ids().collect();
            let data_ids = face_handle.data().attribute_vertices().to_vec();

            // See if there is a vertex id in the face data which is not present in the topology.
            let index = data_ids
                .iter()
                .position(|id| !vert_ids.iter().any(|vid| *vid == *id));
            // If such a vertex exists then it must be the current verex.
            if let Some(i) = index {
                self.mesh.face_data(fid).attribute_vertices_mut()[i] = vid;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::validation::{manifold_state, RedgeManifoldness};
    use crate::wavefront_loader::ObjData;

    use super::*;

    // #[test]
    fn _test_edge_collapse() {
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
        deleter.collapse_edge_and_fix(EdgeId(0));

        let state = manifold_state(deleter.mesh());
        debug_assert!(state == RedgeManifoldness::IsManifold, "{:?}", state);

        let redge = deleter.end_deletion();

        let (vs, fs, _) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/loop_cube.obj");
    }
}
