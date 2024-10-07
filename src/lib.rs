// TODO: For all handles and iterators that can panic, add a fallible API
// wrapper that won't crash.

use std::{collections::BTreeMap, ops::Mul};

pub mod algorithms;
pub mod container_trait;
pub mod edge_handle;
pub mod face_handle;
pub mod hedge_handle;
pub mod helpers;
pub mod iterators;
pub mod mesh_deleter;
pub mod validation;
pub mod vert_handle;

pub use algorithms::quadric_simplification;

use container_trait::{EdgeData, FaceData, PrimitiveContainer, RedgeContainers, VertData};
use edge_handle::EdgeHandle;
use face_handle::FaceHandle;
use hedge_handle::HedgeHandle;
use helpers::{
    _collect_backward_cycle, _collect_forward_cycle, join_radial_cycles, join_vertex_cycles,
    link_face, link_hedges_in_face, remove_edge_from_cycle, split_hedge,
};
use linear_isomorphic::{InnerSpace, RealField};
use validation::{correctness_state, RedgeCorrectness};
use vert_handle::VertHandle;
mod wavefront_loader;

macro_rules! define_id_struct {
    ($name:ident) => {
        #[derive(Clone, Copy, PartialEq, PartialOrd, Ord, Eq, Hash)]
        pub struct $name(usize);
        impl Default for $name {
            fn default() -> Self {
                Self(ABSENT)
            }
        }

        impl $name {
            pub fn is_absent(&self) -> bool {
                self.0 == ABSENT
            }

            pub fn to_index(&self) -> usize {
                self.0
            }
            const ABSENT: Self = Self(ABSENT);
        }

        impl std::fmt::Debug for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                if self.is_absent() {
                    write!(f, "{}(ABSENT)", stringify!($name))
                } else {
                    write!(f, "{}({})", stringify!($name), self.0)
                }
            }
        }
    };
}

pub const ABSENT: usize = usize::MAX;

define_id_struct!(VertId);
define_id_struct!(EdgeId);
define_id_struct!(HedgeId);
define_id_struct!(FaceId);

pub struct Redge<R: RedgeContainers> {
    vert_data: R::VertContainer,
    edge_data: R::EdgeContainer,
    face_data: R::FaceContainer,

    verts_meta: Vec<VertMetaData>,
    edges_meta: Vec<EdgeMetaData>,
    hedges_meta: Vec<HedgeMetaData>,
    faces_meta: Vec<FaceMetaData>,
}

impl<R: RedgeContainers> Redge<R> {
    pub fn new(
        vert_data: R::VertContainer,
        edge_data: R::EdgeContainer,
        face_data: R::FaceContainer,
        faces: impl Iterator<Item = impl Iterator<Item = usize>>,
    ) -> Self {
        // 1. Initialize the setup structures.
        let mut edge_index_map = BTreeMap::<[usize; 2], usize>::new();
        let mut verts_meta = (0..vert_data.len())
            .map(|i| VertMetaData {
                id: VertId(i),
                edge_id: EdgeId(ABSENT),
                is_active: true,
            })
            .collect::<Vec<_>>();
        let mut edges_meta = Vec::new();
        let mut hedges_meta = Vec::new();
        let mut faces_meta = Vec::new();

        // 2. Create all needed edges, faces and hedges. (Set up some of the pointers)
        for face in faces {
            let face_vertices = face.collect::<Vec<_>>();
            let face_index = faces_meta.len();
            faces_meta.push(FaceMetaData {
                id: FaceId(faces_meta.len()),
                is_active: true,
                ..Default::default()
            });
            let active_face = faces_meta.last_mut().unwrap();

            // Remember where we started adding hedges.
            let hedge_cutoff = hedges_meta.len();
            for i in 0..face_vertices.len() {
                let v1 = face_vertices[i];
                let v2 = face_vertices[(i + 1) % face_vertices.len()];
                // Make a canonical representation for the edge by
                // sorting the vertex indices.
                let (cv1, cv2) = if v1 > v2 { (v1, v2) } else { (v2, v1) };

                // Add an edge if it doesn't exist, then fetch a mutable pointer
                // to that edge.
                let edge = if let Some(index) = edge_index_map.get(&[cv1, cv2]) {
                    &mut edges_meta[*index]
                } else {
                    edge_index_map.insert([cv1, cv2], edges_meta.len());
                    edges_meta.push(EdgeMetaData {
                        id: EdgeId(edges_meta.len()),
                        vert_ids: [VertId(cv1), VertId(cv2)],
                        v1_cycle: StarCycleNode {
                            prev_edge: EdgeId(edges_meta.len()),
                            next_edge: EdgeId(edges_meta.len()),
                        },
                        v2_cycle: StarCycleNode {
                            prev_edge: EdgeId(edges_meta.len()),
                            next_edge: EdgeId(edges_meta.len()),
                        },
                        is_active: true,
                        ..Default::default()
                    });
                    edges_meta.last_mut().unwrap()
                };

                // Associate this edge to its endpoints.
                verts_meta[v1].edge_id = edge.id;
                verts_meta[v2].edge_id = edge.id;

                // Add a new hedge.hedge_id
                let hedge_id = HedgeId(hedges_meta.len());
                hedges_meta.push(HedgeMetaData {
                    id: hedge_id,
                    edge_id: edge.id,
                    source_id: VertId(v1),
                    face_id: FaceId(face_index),
                    radial_next_id: hedge_id,
                    radial_prev_id: hedge_id,
                    // The order in which we add hedges to the array is predictable.
                    face_next_id: HedgeId(hedge_cutoff + ((i + 1) % face_vertices.len())),
                    face_prev_id: HedgeId(
                        hedge_cutoff + ((i + face_vertices.len() - 1) % face_vertices.len()),
                    ),
                    is_active: true,
                    ..Default::default()
                });

                active_face.hedge_id = hedge_id;
                edge.hedge_id = hedge_id;
            }
        }

        // 3. Create radial cycles for the hedges.
        for i in 0..hedges_meta.len() {
            let edge = &edges_meta[hedges_meta[i].edge_id.to_index()];
            // This is the active hedge in the edges radial cycle, no need to mutate it directly.
            if edge.hedge_id == hedges_meta[i].id {
                continue;
            }
            // A hedge pointing to itself denotes that it has not been included in its cycle.
            if hedges_meta[i].radial_next_id == hedges_meta[i].radial_prev_id {
                // The prev and next should only ever be equal if they point to
                // the current hedge.
                debug_assert!(
                    hedges_meta[i].radial_next_id == hedges_meta[i].id,
                    "Singleton hedge is incorrect, next_id {:?}, id {:?}",
                    hedges_meta[i].radial_next_id,
                    hedges_meta[i].id
                );

                let edge_id = hedges_meta[i].edge_id;
                debug_assert!(!edge_id.is_absent());

                debug_assert!(!edge_id.is_absent());
                let hedge_id = edges_meta[edge_id.to_index()].hedge_id;
                let next_id = hedges_meta[hedge_id.to_index()].radial_next_id;
                debug_assert!(hedge_id != hedges_meta[i].id);

                // Insert current edge into the radial loop.
                if next_id == hedge_id {
                    // Special case for when there's only one edge in the cycle.
                    hedges_meta[hedge_id.to_index()].radial_next_id = hedges_meta[i].id;
                    hedges_meta[hedge_id.to_index()].radial_prev_id = hedges_meta[i].id;

                    hedges_meta[i].radial_next_id = hedges_meta[hedge_id.to_index()].id;
                    hedges_meta[i].radial_prev_id = hedges_meta[hedge_id.to_index()].id;
                } else {
                    // General logic.
                    hedges_meta[hedge_id.to_index()].radial_next_id = hedges_meta[i].id;
                    hedges_meta[i].radial_prev_id = hedge_id;

                    hedges_meta[next_id.to_index()].radial_prev_id = hedges_meta[i].id;
                    hedges_meta[i].radial_next_id = next_id;
                }
            }
        }

        // 4. For each edge, associate it to its endpoints.
        let mut vertex_edge_associations = vec![Vec::new(); verts_meta.len()];
        for edge in edges_meta.iter() {
            vertex_edge_associations[edge.vert_ids[0].to_index()].push(edge.id);
            vertex_edge_associations[edge.vert_ids[1].to_index()].push(edge.id);
        }

        // 5. Create the edge cycles at the tips of each edge.
        for (vert_index, incident_edges) in vertex_edge_associations.iter().enumerate() {
            for i in 0..incident_edges.len() {
                // Grab the two edges.
                let e1_id = incident_edges[i];
                let e2_id = incident_edges[(i + 1) % incident_edges.len()];

                let vert_id = VertId(vert_index);

                edges_meta[e1_id.to_index()].cycle_mut(vert_id).next_edge = e2_id;
                edges_meta[e2_id.to_index()].cycle_mut(vert_id).prev_edge = e1_id;
            }
        }

        let mesh = Self {
            vert_data,
            edge_data,
            face_data,
            verts_meta,
            edges_meta,
            hedges_meta,
            faces_meta,
        };

        debug_assert!(correctness_state(&mesh) == RedgeCorrectness::Correct);
        mesh
    }

    pub fn vert_count(&self) -> usize {
        self.verts_meta.len()
    }

    pub fn edge_count(&self) -> usize {
        self.edges_meta.len()
    }

    pub fn hedge_count(&self) -> usize {
        self.hedges_meta.len()
    }

    pub fn face_count(&self) -> usize {
        self.faces_meta.len()
    }

    pub fn meta_verts(&self) -> impl Iterator<Item = VertHandle<R>> {
        self.verts_meta.iter().filter_map(|v| {
            if v.is_active {
                Some(self.vert_handle(v.id))
            } else {
                None
            }
        })
    }

    pub fn meta_edges(&self) -> impl Iterator<Item = EdgeHandle<R>> {
        self.edges_meta.iter().filter_map(|e| {
            if e.is_active {
                Some(self.edge_handle(e.id))
            } else {
                None
            }
        })
    }

    pub fn meta_hedges(&self) -> impl Iterator<Item = HedgeHandle<R>> {
        self.hedges_meta.iter().filter_map(|e| {
            if e.is_active {
                Some(self.hedge_handle(e.id))
            } else {
                None
            }
        })
    }

    pub fn meta_faces(&self) -> impl Iterator<Item = FaceHandle<R>> {
        self.faces_meta.iter().filter_map(|f| {
            if f.is_active {
                Some(self.face_handle(f.id))
            } else {
                None
            }
        })
    }

    pub fn vert_handle<'r>(&'r self, id: VertId) -> VertHandle<'r, R> {
        assert!(id.to_index() < self.verts_meta.len());
        VertHandle::new(id, self)
    }

    pub fn edge_handle<'r>(&'r self, id: EdgeId) -> EdgeHandle<'r, R> {
        assert!(id.to_index() < self.edges_meta.len());
        EdgeHandle::new(id, self)
    }

    pub fn hedge_handle<'r>(&'r self, id: HedgeId) -> HedgeHandle<'r, R> {
        assert!(
            id.to_index() < self.hedges_meta.len(),
            "hedge id {:?} is invalid.",
            id
        );
        HedgeHandle::new(id, self)
    }

    pub fn face_handle<'r>(&'r self, id: FaceId) -> FaceHandle<'r, R> {
        assert!(id.to_index() < self.faces_meta.len());
        FaceHandle::new(id, self)
    }

    pub fn vert_data(&mut self, id: VertId) -> &mut VertData<R> {
        debug_assert!(self.verts_meta[id.to_index()].id == id);
        debug_assert!(self.verts_meta[id.to_index()].is_active);
        self.vert_data.get_mut(id.to_index() as u64)
    }

    pub fn edge_data(&mut self, id: EdgeId) -> &mut EdgeData<R> {
        debug_assert!(self.edges_meta[id.to_index()].id == id);
        debug_assert!(self.edges_meta[id.to_index()].is_active);
        self.edge_data.get_mut(id.to_index() as u64)
    }

    pub fn face_data(&mut self, id: FaceId) -> &mut FaceData<R> {
        debug_assert!(self.faces_meta[id.to_index()].id == id);
        debug_assert!(self.faces_meta[id.to_index()].is_active);
        self.face_data.get_mut(id.to_index() as u64)
    }

    /// If you call this in the middle of deletion operations, innactive vertices will still exist.
    /// Since the defragmentation won't be applied until destruction of the deleter, these vertices
    /// will be exported, despite being innactive, to preserve indexing. So if you encounter
    /// ghost or nan vertices in the returned data, check whether those vertices are innactive.
    pub fn to_face_list(&self) -> (Vec<VertData<R>>, Vec<Vec<usize>>) {
        let verts = self.vert_data.iterate().cloned().collect();
        let mut faces = Vec::with_capacity(self.faces_meta.len());

        for face in &self.faces_meta {
            if !face.is_active {
                continue;
            }
            debug_assert!(!face.hedge_id.is_absent());
            let mut current_hedge = &self.hedges_meta[face.hedge_id.to_index()];
            let start_hedge_id = current_hedge.id;

            let mut current_face = Vec::with_capacity(3);
            loop {
                let source_id = current_hedge.source_id;
                current_face.push(source_id.to_index());

                let next_hedge_id = current_hedge.face_next_id;
                debug_assert!(!next_hedge_id.is_absent());

                current_hedge = &self.hedges_meta[next_hedge_id.to_index()];

                if current_hedge.id == start_hedge_id {
                    break;
                }
            }

            faces.push(current_face);
        }

        (verts, faces)
    }

    /// This can only be applied to a non-boundary manifold edge, and only if the incident faces are triangular.
    // If you plan on modifying this function please read `docs/redge.pdf`.
    pub fn flip_edge(&mut self, eid: EdgeId) {
        let edge_handle = self.edge_handle(eid);

        // Get the hedges of the incident faces.
        let hedge = edge_handle.hedge();

        let h1 = hedge.id();
        let h2 = hedge.face_next().id();
        let h3 = hedge.face_prev().id();

        let h4 = hedge.radial_next().id();
        let h5 = hedge.radial_next().face_prev().id();
        let h6 = hedge.radial_next().face_next().id();

        // It's extremely important to access the vertices through the hedge rather than through the edge!
        // This is orientation agnostic the later isn't.
        let v1 = hedge.source().id();
        let v2 = hedge.face_next().source().id();
        let t1 = hedge.face_prev().source().id();
        let t2 = hedge.radial_next().face_prev().source().id();

        // Get the two incident faces before the flip.
        let f1 = hedge.face().id();
        let f2 = hedge.radial_next().face().id();

        let v1_safe_edge = self.vert_handle(v1).pick_different(eid).unwrap().id();
        let v2_safe_edge = self.vert_handle(v2).pick_different(eid).unwrap().id();

        let t1_edge = self.verts_meta[t1.to_index()].edge_id;
        let t2_edge = self.verts_meta[t2.to_index()].edge_id;

        remove_edge_from_cycle(eid, Endpoint::V1, self);
        remove_edge_from_cycle(eid, Endpoint::V2, self);

        // Reconnect hedges in the new configuration and update faces accordingly.
        self.hedges_meta[h1.to_index()].source_id = t1;
        self.hedges_meta[h4.to_index()].source_id = t2;

        link_face(&[h1, h5, h2], f1, self);
        link_face(&[h4, h3, h6], f2, self);

        // Make sure the vertices point to valid edges after the flip.
        self.verts_meta[v1.to_index()].edge_id = v1_safe_edge;
        self.verts_meta[v2.to_index()].edge_id = v2_safe_edge;

        // Update the edge to point to the transversal vertices.
        self.edges_meta[eid.to_index()].vert_ids = [t1, t2];

        join_vertex_cycles(eid, t1_edge, self);
        join_vertex_cycles(eid, t2_edge, self);

        self.verts_meta[t1.to_index()].edge_id = eid;
        self.verts_meta[t2.to_index()].edge_id = eid;
    }

    /// Can only be applied to triangular meshes.
    // If you plan on modifying this function please read `docs/redge.pdf`.
    pub fn split_edge<S>(&mut self, eid: EdgeId) -> VertId
    where
        S: RealField + Mul<VertData<R>, Output = VertData<R>>,
        VertData<R>: InnerSpace<S>,
    {
        let handle = self.edge_handle(eid);

        let mid = (handle.v1().data().clone() + handle.v2().data().clone()) * S::from(0.5).unwrap();
        let vn = self.add_vert(mid);

        let handle = self.edge_handle(eid);
        let hedges_to_split: Vec<_> = handle.hedge().radial_neighbours().map(|h| h.id()).collect();

        let [v1, v2] = handle.vertex_ids();

        let e = eid;
        // TODO: we might need to be smarter about how the data of the new edges is constructed.
        let data = handle.data().clone();
        let en = self.add_edge([v2, vn], data.clone());

        let mut key1 = [v1, vn];
        key1.sort();
        let mut e_hedges = Vec::new();

        let mut key2 = [v2, vn];
        key2.sort();
        let mut en_hedges = Vec::new();

        let mut vn_edges = Vec::new();
        // Split the hedges radially orbitting the edge.
        // Collect the new hedges parallel to the edge, to re-attach them later.
        // Collect the new edges transversal to the edge to attach them later.
        for hid in hedges_to_split {
            let ([h1, h2], et) = split_hedge(vn, hid, self);

            vn_edges.push(et);

            let p1 = self.hedges_meta[h1.to_index()].source_id;
            let p2 = self.hedges_meta[h2.to_index()].source_id;
            let p3 = self.hedge_handle(h2).face_next().source().id();

            let mut h1_key = [p1, p2];
            h1_key.sort();

            // Detect for both hedges which is the one between v_2 and v_n
            // and which is the one between v_n and v_1.
            if h1_key == key1 {
                self.hedges_meta[h1.to_index()].edge_id = e;
                e_hedges.push(h1);
            } else if h1_key == key2 {
                self.hedges_meta[h1.to_index()].edge_id = en;
                en_hedges.push(h1);
            } else {
                panic!();
            }

            let mut h2_key = [p2, p3];
            h2_key.sort();
            if h2_key == key1 {
                self.hedges_meta[h2.to_index()].edge_id = e;
                e_hedges.push(h2);
            } else if h2_key == key2 {
                self.hedges_meta[h2.to_index()].edge_id = en;
                en_hedges.push(h2);
            } else {
                panic!();
            }
        }

        self.edges_meta[e.to_index()].hedge_id = e_hedges[0];
        self.edges_meta[en.to_index()].hedge_id = en_hedges[0];
        vn_edges.push(e);
        vn_edges.push(en);

        let safe_edge = self
            .edge_handle(eid)
            .v2()
            .star_edges()
            .find(|e| e.id() != eid)
            .expect("Trying to split an edge with degenerate topology.")
            .id();

        // Reset the vertex cycle at the opposite endpoint of the current edge.
        remove_edge_from_cycle(eid, Endpoint::V2, self);
        self.verts_meta[v2.to_index()].edge_id = safe_edge;
        *self.edges_meta[e.to_index()].at(v2) = vn;

        let cycle = self.edges_meta[e.to_index()].cycle_mut(vn);
        cycle.next_edge = e;
        cycle.prev_edge = e;

        // Join the radial cycles of the newly created edges after the split.
        for i in 0..e_hedges.len() {
            let h1 = e_hedges[i];
            let h2 = e_hedges[(i + 1) % e_hedges.len()];

            self.hedges_meta[h1.to_index()].radial_next_id = h2;
            self.hedges_meta[h2.to_index()].radial_prev_id = h1;
        }

        for i in 0..en_hedges.len() {
            let h1 = en_hedges[i];
            let h2 = en_hedges[(i + 1) % en_hedges.len()];

            self.hedges_meta[h1.to_index()].radial_next_id = h2;
            self.hedges_meta[h2.to_index()].radial_prev_id = h1;
        }

        // Create a new cycle of edges around the new vertex.
        for i in 0..vn_edges.len() - 1 {
            let e1 = vn_edges[i];
            let e2 = vn_edges[i + 1];

            join_vertex_cycles(e1, e2, self);
        }

        self.verts_meta[vn.to_index()].edge_id = en;

        let v2_edge = self.verts_meta[v2.to_index()].edge_id;
        join_vertex_cycles(v2_edge, en, self);

        debug_assert!(correctness_state(&self) == RedgeCorrectness::Correct);

        vn
    }

    pub(crate) fn add_vert(&mut self, data: VertData<R>) -> VertId {
        self.vert_data.push(data);
        let id = VertId(self.verts_meta.len() as usize);
        self.verts_meta.push(VertMetaData {
            id,
            edge_id: EdgeId::ABSENT,
            is_active: true,
        });

        id
    }

    pub(crate) fn add_edge(&mut self, vert_ids: [VertId; 2], data: EdgeData<R>) -> EdgeId {
        let id = EdgeId(self.edges_meta.len() as usize);
        self.edge_data.push(data);
        self.edges_meta.push(EdgeMetaData {
            id,
            vert_ids,
            hedge_id: HedgeId::ABSENT,
            v1_cycle: StarCycleNode {
                next_edge: id,
                prev_edge: id,
            },
            v2_cycle: StarCycleNode {
                next_edge: id,
                prev_edge: id,
            },
            is_active: true,
        });

        id
    }

    /// Warning: Unlike `add_vert` and `add_edge`, this function will leave
    /// the redge in an invalid state because the hedge and edge pointers won't be initialized.
    /// Use carefully.
    pub(crate) fn add_hedge(&mut self, source_id: VertId) -> HedgeId {
        let id = HedgeId(self.hedges_meta.len() as usize);
        self.hedges_meta.push(HedgeMetaData {
            id,
            edge_id: EdgeId::ABSENT,
            radial_next_id: id,
            radial_prev_id: id,
            face_next_id: HedgeId::ABSENT,
            face_prev_id: HedgeId::ABSENT,
            source_id,
            face_id: FaceId::ABSENT,
            is_active: true,
        });

        id
    }

    /// Warning: Unlike `add_vert` and `add_edge`, this function will leave
    /// the redge in an invalid state because the face pointers won't point to anything. Use very carefully.
    pub(crate) fn add_face(&mut self, data: FaceData<R>) -> FaceId {
        self.face_data.push(data);
        let id = FaceId(self.faces_meta.len() as usize);

        self.faces_meta.push(FaceMetaData {
            id,
            hedge_id: HedgeId::ABSENT,
            is_active: true,
        });

        id
    }
}

#[derive(Default, Debug, Clone)]
struct VertMetaData {
    id: VertId,
    /// Points to any edge that touches this vertex.
    edge_id: EdgeId,
    /// Whether this element is being used (set to false on removal).
    is_active: bool,
}
#[derive(Default, Debug, Clone)]
struct EdgeMetaData {
    id: EdgeId,
    /// The two endpoints of the edge. It is an invariant that
    /// they are sorted by their id.
    vert_ids: [VertId; 2],
    /// Any hedge on a face alongside this edge.
    hedge_id: HedgeId,
    /// Edges touching the v1 endpoint.
    v1_cycle: StarCycleNode,
    /// Edges touching the v2 endpoint.
    v2_cycle: StarCycleNode,
    /// Whether this element is being used (set to false on removal).
    is_active: bool,
}

impl EdgeMetaData {
    pub(crate) fn cycle(&self, vert_id: VertId) -> &StarCycleNode {
        debug_assert!(self.is_active);
        if vert_id == self.vert_ids[0] {
            return &self.v1_cycle;
        } else if vert_id == self.vert_ids[1] {
            return &self.v2_cycle;
        } else {
            panic!(
                "Requested vert {:?} in edge ({:?}, {:?})",
                vert_id, self.vert_ids[0], self.vert_ids[1],
            )
        }
    }

    pub(crate) fn cycle_mut(&mut self, vert_id: VertId) -> &mut StarCycleNode {
        if vert_id == self.vert_ids[0] {
            return &mut self.v1_cycle;
        } else if vert_id == self.vert_ids[1] {
            return &mut self.v2_cycle;
        } else {
            panic!("Vertex {:?} does not exist in edge {:?}.", vert_id, self.id)
        }
    }

    pub(crate) fn at(&mut self, vert_id: VertId) -> &mut VertId {
        if vert_id == self.vert_ids[0] {
            return &mut self.vert_ids[0];
        } else if vert_id == self.vert_ids[1] {
            return &mut self.vert_ids[1];
        } else {
            panic!("{:?}", vert_id)
        }
    }

    pub(crate) fn topology_intersection(&self, other: &EdgeMetaData) -> Option<VertId> {
        let vids = self.vert_ids;
        let oids = other.vert_ids;

        if vids[0] == oids[0] || vids[0] == oids[1] {
            return Some(vids[0]);
        } else if vids[1] == oids[0] || vids[1] == oids[1] {
            return Some(vids[1]);
        }

        return None;
    }
}

/// This is the radial part of the edge. It's part of the cycle of all the
/// directed edges of all faces that touch an edge. The `radial_*` pointers
/// refer to directed edges that are parallel to each other. For example,
/// if 3 faces meet at the edge, then the pointers form a doubly linked
/// list of 3 elements that radially orbits the edge.
/// The `face_*` pointers are the cycle of directed edges along a face, obeying
/// its orientation.
#[derive(Default, Debug, Clone)]
struct HedgeMetaData {
    /// Unique identfier for the hedge.
    id: HedgeId,
    /// Vertex from which this hedge starts at.
    source_id: VertId,
    /// Edge that is parallel to this hedge and which is contained by the face
    /// of this hedge.
    edge_id: EdgeId,
    /// Next hedge in the parallel cycle around the edge.
    radial_next_id: HedgeId,
    /// Prev hedge in the parallel cycle around the edge.
    radial_prev_id: HedgeId,
    /// Next hedge in the face loop.
    face_next_id: HedgeId,
    /// Prev hedge in the face loop.
    face_prev_id: HedgeId,
    /// The face that contains this hedge.
    face_id: FaceId,
    /// Whether this element is being used (set to false on removal).
    is_active: bool,
}

#[derive(Default, Debug, Clone)]
struct FaceMetaData {
    id: FaceId,
    hedge_id: HedgeId,
    /// Whether this element is being used (set to false on removal).
    is_active: bool,
}

#[derive(Default, Debug, Clone)]
pub struct StarCycleNode {
    prev_edge: EdgeId,
    next_edge: EdgeId,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Endpoint {
    V1,
    V2,
}

#[cfg(test)]
mod tests {
    use nalgebra::Vector3;
    use validation::{manifold_state, RedgeManifoldness};
    use wavefront_loader::ObjData;

    use super::*;

    #[test]
    fn test_edge_flip() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/loop_cube.obj");

        let mut redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/before_flip.obj");

        redge.flip_edge(EdgeId(0));

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/after_flip.obj");
    }

    #[test]
    fn test_edge_split() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/loop_cube.obj");

        // Convert to a type that admits arithmetic transformations.
        let vertices: Vec<_> = vertices
            .iter()
            .map(|v| Vector3::new(v[0], v[1], v[2]))
            .collect();

        let mut redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/before_flip.obj");

        redge.split_edge(EdgeId(0));

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/after_flip.obj");
    }

    #[test]
    fn test_redge_init() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/tetrahedron.obj");

        let redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        debug_assert!(manifold_state(&redge) == RedgeManifoldness::IsManifold,);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/tetrahedron.obj");

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

        debug_assert!(manifold_state(&redge) == RedgeManifoldness::IsManifold);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/loop_cube.obj");

        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/non_manifold_tet.obj");

        let redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let state = manifold_state(&redge);
        debug_assert!(state != RedgeManifoldness::IsManifold, "{:?}", state);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/non_manifold_tet.obj");

        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/triangle.obj");

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

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/triangle.obj");

        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/triforce.obj");

        let redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let state = manifold_state(&redge);
        debug_assert!(state != RedgeManifoldness::IsManifold, "{:?}", state);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/triforce.obj");

        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/triple_triangle.obj");

        let redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let state = manifold_state(&redge);
        debug_assert!(state != RedgeManifoldness::IsManifold, "{:?}", state);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/triple_triangle.obj");
    }
}
