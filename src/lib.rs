// TODO: For all handles and iterators that can panic, add a fallible API
// warapper that won't crash.

use std::collections::BTreeMap;

pub mod container_trait;
pub mod edge_handle;
pub mod face_handle;
pub mod hedge_handle;
pub mod iterators;
pub mod validation;
pub mod vert_handle;
use container_trait::PrimitiveContainer;
use edge_handle::EdgeHandle;
use face_handle::FaceHandle;
use hedge_handle::HedgeHandle;
use vert_handle::VertHandle;
mod wavefront_loader;

macro_rules! define_id_struct {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Ord, Eq)]
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
        }
    };
}

pub const ABSENT: usize = usize::MAX;

define_id_struct!(VertId);
define_id_struct!(EdgeId);
define_id_struct!(HedgeId);
define_id_struct!(FaceId);

pub struct Redge<VContainer, EContainer, FContainer>
where
    VContainer: PrimitiveContainer,
    EContainer: PrimitiveContainer,
    FContainer: PrimitiveContainer,
{
    vert_data: VContainer,
    edge_data: EContainer,
    face_data: FContainer,

    verts_meta: Vec<VertMetaData>,
    edges_meta: Vec<EdgeMetaData>,
    hedges_meta: Vec<HedgeMetaData>,
    faces_meta: Vec<FaceMetaData>,
}

impl<VContainer, EContainer, FContainer> Redge<VContainer, EContainer, FContainer>
where
    VContainer: PrimitiveContainer,
    EContainer: PrimitiveContainer,
    FContainer: PrimitiveContainer,
{
    pub fn new(
        vert_data: VContainer,
        edge_data: EContainer,
        face_data: FContainer,
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
            }
        }

        // 3. Enforce that there are no singleton radial hedges.
        for i in 0..hedges_meta.len() {
            // If a hedge has no other elements in its cycle, we will make a new hedge
            // in the opposite direction and add it to the cycle.
            if hedges_meta[i].radial_next_id == hedges_meta[i].radial_prev_id {
                // The prev and next should only ever be equal if they point to
                // the same hedge.
                debug_assert!(
                    hedges_meta[i].radial_next_id == hedges_meta[i].id,
                    "Singleton hedge is incorrect, next_id {:?}, id {:?}",
                    hedges_meta[i].radial_next_id,
                    hedges_meta[i].id
                );

                let edge_id = hedges_meta[i].edge_id;
                debug_assert!(!edge_id.is_absent());
                let source_id = hedges_meta[i].source_id;
                debug_assert!(!source_id.is_absent());

                // The source of the pair is the destination of the current
                // hedge.
                let [v1, v2] = edges_meta[edge_id.to_index()].vert_ids;
                let dest = if source_id == v1 { v2 } else { v1 };

                let hedge_id = HedgeId(hedges_meta.len());
                hedges_meta.push(HedgeMetaData {
                    id: hedge_id,
                    source_id: dest,
                    is_active: true,
                    ..Default::default()
                });

                attach_hedge_to_edge(
                    &mut edges_meta[edge_id.to_index()],
                    hedge_id,
                    &mut hedges_meta,
                );
            }
        }

        // 4. For each edge, associate it to its endpoints.
        let mut vertex_edge_associations = vec![Vec::new(); verts_meta.len()];
        for edge in edges_meta.iter() {
            vertex_edge_associations[edge.vert_ids[0].to_index()].push((Endpoint::V1, edge.id));
            vertex_edge_associations[edge.vert_ids[1].to_index()].push((Endpoint::V2, edge.id));
        }

        // 5. Create the edge cycles at the tips of each edge.
        for (vert_index, incident_edges) in vertex_edge_associations.iter().enumerate() {
            let vert_id = VertId(vert_index);
            for i in 0..incident_edges.len() {
                // Grab the two edges.
                let (endpoint1, e1_id) = incident_edges[i];
                let (endpoint2, e2_id) = incident_edges[(i + 1) % incident_edges.len()];
                let (e1, e2) = get_disjoint(&mut edges_meta, e1_id.to_index(), e2_id.to_index());

                // Get a handle to the correct cycle for each edge.
                let cycle1 = match endpoint1 {
                    Endpoint::V1 => &mut e1.v1_cycle,
                    Endpoint::V2 => &mut e1.v2_cycle,
                };
                let cycle2 = match endpoint2 {
                    Endpoint::V1 => &mut e2.v1_cycle,
                    Endpoint::V2 => &mut e2.v2_cycle,
                };

                // Link the edge endpoint chain in the given order.
                cycle1.next_edge = e2.id;
                cycle2.prev_edge = e1.id;
            }
        }

        Self {
            vert_data,
            edge_data,
            face_data,
            verts_meta,
            edges_meta,
            hedges_meta,
            faces_meta,
        }
    }

    pub fn vert_handle<'r>(
        &'r self,
        id: VertId,
    ) -> VertHandle<'r, VContainer, EContainer, FContainer> {
        assert!(id.to_index() < self.verts_meta.len());
        VertHandle::new(id, self)
    }

    pub fn edge_handle<'r>(
        &'r self,
        id: EdgeId,
    ) -> EdgeHandle<'r, VContainer, EContainer, FContainer> {
        assert!(id.to_index() < self.edges_meta.len());
        EdgeHandle::new(id, self)
    }

    pub fn hedge_handle<'r>(
        &'r self,
        id: HedgeId,
    ) -> HedgeHandle<'r, VContainer, EContainer, FContainer> {
        assert!(id.to_index() < self.hedges_meta.len());
        HedgeHandle::new(id, self)
    }

    pub fn face_handle<'r>(
        &'r self,
        id: FaceId,
    ) -> FaceHandle<'r, VContainer, EContainer, FContainer> {
        assert!(id.to_index() < self.faces_meta.len());
        FaceHandle::new(id, self)
    }

    pub fn to_face_list(&self) -> (Vec<VContainer::PrimitiveData>, Vec<Vec<usize>>) {
        let verts = self.vert_data.iterate().cloned().collect();
        let mut faces = Vec::with_capacity(self.faces_meta.len());

        for face in &self.faces_meta {
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
    /// The two endpoints of the edge.
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

/// This is the radial part of the edge. It's part of the cycle of all the
/// directed edges of all faces that touch an edge. The `radial_*` pointers
/// refer to directed edges that are parallel to each other. For example,
/// if 3 faces meet at the edge, then the pointers form a doubly linked
/// list of 3 elements that radially orbits the edge.
/// The `face_*` pointers are the cycle of directed edges along a face, obeying
/// its orientation.
#[derive(Default, Debug, Clone)]
struct HedgeMetaData {
    id: HedgeId,
    source_id: VertId,
    edge_id: EdgeId,
    /// Next hedge in the parallel cycle around the edge.
    radial_next_id: HedgeId,
    /// Prev hedge in the parallel cycle around the edge.
    radial_prev_id: HedgeId,
    /// Next hedge in the face.
    face_next_id: HedgeId,
    /// Prev hedge in the face.
    face_prev_id: HedgeId,
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

#[derive(Debug, Clone, Copy)]
enum Endpoint {
    V1,
    V2,
}

#[inline]
fn get_disjoint<'a, T>(x: &'a mut [T], a: usize, b: usize) -> (&'a mut T, &'a mut T) {
    assert_ne!(a, b);
    assert!(a < x.len());
    assert!(b < x.len());
    let ptr = x.as_mut_ptr();
    (unsafe { &mut *ptr.add(a) }, unsafe { &mut *ptr.add(b) })
}

#[inline]
fn attach_hedge_to_edge(
    edge_meta: &mut EdgeMetaData,
    new_hedge_id: HedgeId,
    hedges: &mut Vec<HedgeMetaData>,
) {
    if edge_meta.hedge_id.is_absent() {
        edge_meta.hedge_id = new_hedge_id;
        hedges[new_hedge_id.to_index()].radial_next_id = new_hedge_id;
        hedges[new_hedge_id.to_index()].radial_prev_id = new_hedge_id;
    } else {
        insert_hedge_in_cycle(edge_meta.hedge_id, new_hedge_id, hedges);
    }
}
/// Insert `new` into the doubly linked list containing `old_id` between `old_id` and its next.
#[inline]
fn insert_hedge_in_cycle(old_id: HedgeId, new_id: HedgeId, hedges: &mut Vec<HedgeMetaData>) {
    let old_next_id = hedges[old_id.to_index()].radial_next_id;

    link_hedge_in_face(old_id, new_id, hedges);
    link_hedge_in_face(new_id, old_next_id, hedges);
}

#[inline]
fn link_hedge_in_face(prev_id: HedgeId, next_id: HedgeId, hedges: &mut Vec<HedgeMetaData>) {
    hedges[prev_id.to_index()].face_next_id = next_id;
    hedges[next_id.to_index()].face_prev_id = prev_id;
}

#[cfg(test)]
mod tests {
    use validation::{manfiold_state, RedgeManifoldness};
    use wavefront_loader::ObjData;

    use super::*;

    #[test]
    fn test_redge_init() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/tetrahedron.obj");

        let redge = Redge::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let manifold_state = manfiold_state(&redge);
        debug_assert!(
            manifold_state == RedgeManifoldness::IsManifold,
            "{:?}",
            manifold_state
        );

        let (vs, fs) = redge.to_face_list();

        ObjData::export(&(&vs, &fs), "out/tetrahedron.obj");

        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/hedge_cube.obj");

        let redge = Redge::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let (vs, fs) = redge.to_face_list();

        ObjData::export(&(&vs, &fs), "out/hedge_cube.obj");
    }
}
