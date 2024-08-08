use std::collections::BTreeSet;

use crate::{
    container_trait::RedgeContainers,
    face_handle,
    helpers::{
        disable_edge_meta, disable_face_meta, disable_hedge_meta, disable_vert_meta,
        join_radial_cycles, join_vertex_cycles, remove_edge_from_cycle, remove_hedge_from_face,
        remove_hedge_from_radial,
    },
    EdgeId, FaceId, HedgeId, Redge, StarCycleNode, VertId,
};

pub struct MeshDeleter<R: RedgeContainers> {
    mesh: Redge<R>,
}

impl<R: RedgeContainers> MeshDeleter<R> {
    pub fn start_deletion(mesh: Redge<R>) -> Self {
        Self { mesh }
    }

    pub fn end_deletion(self) -> Redge<R> {
        self.mesh
    }

    pub fn collapse_edge(&mut self) {}

    /// Collapse a face's hedge during an edge collapse.
    fn edge_face_remove(&mut self, hedge_id: HedgeId) {
        let hedge_handle = self.mesh.hedge_handle(hedge_id);
        let f = hedge_handle.face().id();

        let hn = hedge_handle.face_next().id();
        let hp = hedge_handle.face_prev().id();

        let h0 = hedge_handle.id();

        let v1 = hedge_handle.face_next().source().id();

        remove_hedge_from_face(h0, &mut self.mesh);
        remove_hedge_from_radial(h0, &mut self.mesh);

        self.mesh.faces_meta[f.to_index()].hedge_id = hp;

        disable_hedge_meta(h0, &mut self.mesh);

        self.mesh.hedges_meta[hn.to_index()].source_id = v1;
    }

    /// Collapse an edge (just the edge) during an edge collapse.
    fn edge_edge_remove(&mut self, edge_id: EdgeId) {
        let edge_handle = self.mesh.edge_handle(edge_id);
        let v1 = edge_handle.v1().id();
        let v2 = edge_handle.v2().id();
        let c1 = edge_handle.vertex_cycle_pointers_v1();
        let c2 = edge_handle.vertex_cycle_pointers_v2();
        let e0 = edge_handle.id();

        // Select a safe edge that will survive the collapse.
        let mut iterator = edge_handle.v1().iter_star_edges();
        let e_safe = match iterator.next() {
            None => EdgeId::new_absent(),
            Some(_) => iterator.next().unwrap().id(),
        };

        let v2_neighbours: Vec<_> = self
            .mesh
            .vert_handle(v2)
            .iter_star_edges()
            .map(|v| v.id())
            .collect();

        remove_edge_from_cycle(e0, crate::Endpoint::V1, &mut self.mesh);
        remove_edge_from_cycle(e0, crate::Endpoint::V2, &mut self.mesh);

        // Make all old neighbours of v2 now point to v1.
        for id in v2_neighbours {
            let index = self.mesh.edges_meta[id.to_index()].local_index(v2);
            self.mesh.edges_meta[id.to_index()].vert_ids[index] = v1;
        }

        join_vertex_cycles(c1, c2, v1, &mut self.mesh);

        disable_edge_meta(e0, &mut self.mesh);
        disable_vert_meta(v2, &mut self.mesh);

        self.mesh.verts_meta[v1.to_index()].edge_id = e_safe;
    }

    /// This assumes the face is degenerate and doesn't check
    /// be careful when calling it.
    fn remove_degenerate_face(&mut self, face_id: FaceId) {
        let face = self.mesh.face_handle(face_id);
        let f = face.id();

        let h1 = face.hedge().face_prev().id();
        let h2 = face.hedge().face_next().id();

        let v1 = face.hedge().edge().v1().id();
        let v2 = face.hedge().edge().v2().id();

        let e1 = face.hedge().edge().id();
        let e2 = face.hedge().face_next().edge().id();

        let h1_r = face.hedge().radial_next().id();
        let h2_r = face.hedge().face_next().radial_next().id();

        // Collect all half edges around e2 that are neither
        // h1 nor h2.
        let surviving: BTreeSet<_> = self
            .mesh
            .edge_handle(e2)
            .hedge()
            .radial_neighbours()
            .filter_map(|e| {
                if e.id() != h2 && e.id() != h1 {
                    Some(e.id())
                } else {
                    None
                }
            })
            .collect();

        remove_hedge_from_radial(h1, &mut self.mesh);
        remove_hedge_from_radial(h2, &mut self.mesh);

        remove_edge_from_cycle(e2, crate::Endpoint::V1, &mut self.mesh);
        remove_edge_from_cycle(e2, crate::Endpoint::V2, &mut self.mesh);

        join_radial_cycles(h1_r, h2_r, &mut self.mesh);

        for survivor in surviving {
            self.mesh.hedges_meta[survivor.to_index()].edge_id = e1;
        }

        disable_face_meta(f, &mut self.mesh);
        disable_edge_meta(e2, &mut self.mesh);

        self.mesh.edges_meta[e1.to_index()].hedge_id = h1_r;

        self.mesh.verts_meta[v1.to_index()].edge_id = e1;
        self.mesh.verts_meta[v2.to_index()].edge_id = e1;
    }
}
