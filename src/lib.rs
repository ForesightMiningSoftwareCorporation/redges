use std::collections::BTreeMap;

use linear_isomorphic::*;

pub mod container_trait;
use container_trait::PrimitiveContainer;
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
define_id_struct!(LoopId);
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
    loops_meta: Vec<LoopMetaData>,
    faces_meta: Vec<FaceMetaData>,
}

impl<VContainer, EContainer, FContainer> Redge<VContainer, EContainer, FContainer>
where
    VContainer: PrimitiveContainer,
    EContainer: PrimitiveContainer,
    FContainer: PrimitiveContainer,
{
    pub fn new<I>(
        vert_data: VContainer,
        edge_data: EContainer,
        face_data: FContainer,
        faces: impl Iterator<Item = impl Iterator<Item = usize>>,
    ) {
        // 1. Initialize the setup structures.
        let mut edge_index_map = BTreeMap::<[usize; 2], usize>::new();
        let mut verts_meta = (0..vert_data.len())
            .map(|i| VertMetaData {
                id: VertId(i),
                edge_id: EdgeId(ABSENT),
            })
            .collect::<Vec<_>>();
        let mut edges_meta = Vec::new();
        let mut loops_meta = Vec::new();
        let mut faces_meta = Vec::new();

        // 2. Create all needed edges, faces and loops. (Set up some of the pointers)
        for face in faces {
            let face_vertices = face.collect::<Vec<_>>();
            let face_index = faces_meta.len();
            faces_meta.push(FaceMetaData {
                id: FaceId(faces_meta.len()),
                ..Default::default()
            });
            let active_face = faces_meta.last_mut().unwrap();

            // Remember where we started adding loops.
            let loop_cutoff = loops_meta.len();
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
                        vert_ids: [VertId(cv1), VertId(cv2)],
                        ..Default::default()
                    });
                    edges_meta.last_mut().unwrap()
                };

                // Associate this edge to its endpoints.
                verts_meta[v1].edge_id = edge.id;
                verts_meta[v2].edge_id = edge.id;

                // Add a new loop.
                let loop_id = LoopId(loops_meta.len());
                loops_meta.push(LoopMetaData {
                    id: loop_id,
                    edge_id: edge.id,
                    source_id: VertId(v1),
                    face_id: FaceId(face_index),
                    ..Default::default()
                });
                attach_loop_to_edge(edge, loop_id, &mut loops_meta);

                active_face.loop_id = loop_id;
            }

            // Create a doubly linked list for the loops around the face, in the
            // order in which the edges were declared.
            for i in 0..face_vertices.len() {
                let prev_id = LoopId(loop_cutoff + i);
                let next_id = LoopId(loop_cutoff + (i % face_vertices.len()));

                link_loop_in_face(prev_id, next_id, &mut loops_meta);
            }
        }

        // 3. Enforce that there are no singleton radial loops.
        for i in 0..loops_meta.len() {
            // If a loop has no other elements in its cycle, we will make a new loop
            // in the opposite direction and add it to the cycle.
            if loops_meta[i].radial_next_id == loops_meta[i].radial_prev_id {
                // The prev and next should only ever be equal if they point to
                // the same loop.
                debug_assert!(loops_meta[i].radial_next_id == loops_meta[i].id);

                let edge_id = loops_meta[i].edge_id;
                debug_assert!(!edge_id.is_absent());
                let source_id = loops_meta[i].source_id;
                debug_assert!(!source_id.is_absent());

                // The source of the pair is the destination of the current
                // loop.
                let [v1, v2] = edges_meta[edge_id.to_index()].vert_ids;
                let dest = if source_id == v1 { v2 } else { v1 };

                let loop_id = LoopId(loops_meta.len());
                loops_meta.push(LoopMetaData {
                    id: loop_id,
                    source_id: dest,
                    ..Default::default()
                });

                attach_loop_to_edge(
                    &mut edges_meta[edge_id.to_index()],
                    loop_id,
                    &mut loops_meta,
                );
                link_loop_in_face(LoopId(i), loop_id, &mut loops_meta);
            }
        }
    }
}

#[inline]
fn attach_loop_to_edge(
    edge_meta: &mut EdgeMetaData,
    new_loop_id: LoopId,
    loops: &mut Vec<LoopMetaData>,
) {
    if edge_meta.loop_id.is_absent() {
        edge_meta.loop_id = new_loop_id;
        loops[new_loop_id.to_index()].radial_next_id = new_loop_id;
        loops[new_loop_id.to_index()].radial_prev_id = new_loop_id;
    } else {
        insert_loop_in_cycle(edge_meta.loop_id, new_loop_id, loops);
    }
}
/// Insert `new` into the doubly linked list containing `old_id` between `old_id` and its next.
#[inline]
fn insert_loop_in_cycle(old_id: LoopId, new_id: LoopId, loops: &mut Vec<LoopMetaData>) {
    let old_next_id = loops[old_id.to_index()].radial_next_id;

    link_loop_in_face(old_id, new_id, loops);
    link_loop_in_face(new_id, old_next_id, loops);
}

#[inline]
fn link_loop_in_face(prev_id: LoopId, next_id: LoopId, loops: &mut Vec<LoopMetaData>) {
    loops[prev_id.to_index()].radial_next_id = next_id;
    loops[next_id.to_index()].radial_prev_id = prev_id;
}

#[derive(Default, Debug, Clone)]
struct VertMetaData {
    id: VertId,
    /// Points to any edge that touches this vertex.
    edge_id: EdgeId,
}
#[derive(Default, Debug, Clone)]
struct EdgeMetaData {
    id: EdgeId,
    /// The two endpoints of the edge.
    vert_ids: [VertId; 2],
    /// Any loop on a face alongside this edge.
    loop_id: LoopId,
    /// Edges touching the v1 endpoint.
    v1_edges: StarCycleNode,
    /// Edges touching the v2 endpoint.
    v2_edges: StarCycleNode,
}

/// This is the radial part of the edge. It's part of the cycle of all the
/// directed edges of all faces that touch an edge. The `radial_*` pointers
/// refer to directed edges that are parallel to each other. For example,
/// if 3 faces meet at the edge, then the pointers form a doubly linked
/// list of 3 elements that radially orbites the edge.
/// The `face_*` pointers are the cycle of directed edges along a face, obeying
/// its orientation.
#[derive(Default, Debug, Clone)]
struct LoopMetaData {
    id: LoopId,
    source_id: VertId,
    edge_id: EdgeId,
    /// Next loop in the parallel cycle around the edge.
    radial_next_id: LoopId,
    /// Prev loop in the parallel cycle around the edge.
    radial_prev_id: LoopId,
    /// Next loop in the face.
    face_next_id: LoopId,
    /// Prev loop in the face.
    face_prev_id: LoopId,
    face_id: FaceId,
}

#[derive(Default, Debug, Clone)]
struct FaceMetaData {
    id: FaceId,
    loop_id: LoopId,
}

#[derive(Default, Debug, Clone)]
struct StarCycleNode {
    prev_edege: EdgeId,
    next_edge: EdgeId,
}
