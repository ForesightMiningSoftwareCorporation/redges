//! Functions to verify topological and geometrical invariants and validity.

// The functions defined here should remain functions, don't refactor them to methods.
// The reason is that, we will likely want to implement trait interfaces
// to abstract meshes. These methods should work for those abstractions.
// If we couple them to the objects it will be harder to do this.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Debug,
    marker::PhantomData,
};

use linear_isomorphic::{InnerSpace, RealField};
use num::{Bounded, Signed};
use rstar::RTree;

use crate::{
    container_trait::{PrimitiveContainer, RedgeContainers, VertData},
    helpers::check_edge_vertex_cycles,
    EdgeId, Endpoint, HedgeId, Redge, VertId,
};

/// Marks the relationship between a bad element and the element that was inspected.
#[derive(Debug, PartialEq, Eq)]
pub enum TopologicalRelation {
    /// This element is next to the inspected element.
    Next,
    /// This element is previous to the inspected element.
    Prior,
}

// TODO: For each function in this file there will be missing checks, add them as
// needed until they are entirely trustworthy.
/// Error reporting values for the first degeneracy found in a radial edge.
#[derive(Debug, PartialEq, Eq)]
pub enum RedgeCorrectness {
    /// Radial edge is adequate.
    Correct,
    /// Counts between metadata and data do not match.
    MismatchingArrayLengths,
    /// An active element points to an invalid vertex at the specified internal index.
    InvalidVert(usize),
    /// An active element points to an invalid edge at the specified internal index.
    InvalidEdge(usize, EdgeCorrectness),
    /// An active half edge points to an invalid vertex at the specified internal index.
    InvalidHedge(usize, HedgeCorrectness),
    /// An active element points to an invalid face at the specified internal index.
    InvalidFace(usize),
    /// The vertex cycles of two different edges incident at the same common vertex do not match.
    VertexCyclesDontMatch(VertId, BTreeSet<EdgeId>, BTreeSet<EdgeId>),
}

/// Error reporting values for the first degeneracy found in an edge.
#[derive(Debug, PartialEq, Eq)]
pub enum EdgeCorrectness {
    /// Edge is adequate.
    Correct,
    /// An active element is pointing to an absent edge.
    EdgeIsAbsent,
    /// The id of this edge does not match its internal index.
    IdAndIndexMismatch,
    /// One of the endpoints is pointing to an absent vertex.
    AbsentEndpoint(Endpoint),
    /// One of the radial cycles is empty (this should not happen even if there is a single isolated
    /// edge in the entire mesh).
    CyclePointsToAbsent(Endpoint),
    /// One of the endpoint cycles is not properly set up.
    CycleIsBroken(Endpoint, TopologicalRelation, VertId),
    /// Both endpoints of this edge point to the same underlying vertex.
    EdgeWithOnlyOnePoint,
    /// The hedge pointer points to an invalid/removed half edge.
    InvalidHedgePointer(HedgeId),
    /// The half edge this edge points to, points to a different edge.
    HedgePointsToDifferentEdge,
}

/// Error reporting values for the first degeneracy found in a half edge.
#[derive(Debug, PartialEq, Eq)]
pub enum HedgeCorrectness {
    /// Half edge is adequate.
    Correct,
    /// The id of this half edge does not correspond to its internal index.
    IdAndIndexMismatch,
    /// The loop this half edge points to dos not form a simple cycle around the face.
    FaceLoopChainIsBroken,
    /// The loop this half edge points to dos not form a simple cycle orbiting the edge.
    RadialChainIsBroken,
    /// The source vertex id points to an absent/invalid vertex.
    SourceIsAbsent,
    /// The edge if points to an absent/invalid vertex.
    EdgeIsAbsent,
    /// The radial chain has more than 100 elements. This is not, strictly speaking, an error
    /// but it's extremely suspicious.
    LongRadialChain(HedgeId),
}

/// If this returns `Correct` then it is safe to create handles.
// If the above condition is found to be false ping Makogan and tell him
// to fix asap.
pub fn correctness_state<R: RedgeContainers>(mesh: &Redge<R>) -> RedgeCorrectness {
    if mesh.vert_data.len() != mesh.verts_meta.len()
        || mesh.edge_data.len() != mesh.edges_meta.len() && mesh.edge_data.len() != 0
        || mesh.face_data.len() != mesh.faces_meta.len() && mesh.face_data.len() != 0
    {
        return RedgeCorrectness::MismatchingArrayLengths;
    }

    // Validate verts.
    for (i, vert) in mesh.verts_meta.iter().enumerate() {
        if !vert.is_active {
            continue;
        }
        if vert.id.to_index() != i {
            return RedgeCorrectness::InvalidVert(i);
        }

        if vert.edge_id.is_absent() {
            continue;
        }

        if vert.edge_id.to_index() >= mesh.edges_meta.len() {
            return RedgeCorrectness::InvalidVert(i);
        }

        let edge_handle = mesh.edge_handle(vert.edge_id);

        let [v1, v2] = edge_handle.vertex_ids();
        if v1 != vert.id && v2 != vert.id {
            return RedgeCorrectness::InvalidVert(i);
        }
    }

    // Validate edges.
    for (i, edge) in mesh.edges_meta.iter().enumerate() {
        if !edge.is_active {
            continue;
        }
        if edge.id.to_index() != i {
            return RedgeCorrectness::InvalidEdge(i, EdgeCorrectness::IdAndIndexMismatch);
        }
        if edge.vert_ids[0].is_absent() {
            return RedgeCorrectness::InvalidEdge(i, EdgeCorrectness::AbsentEndpoint(Endpoint::V1));
        }
        if edge.vert_ids[1].is_absent() {
            return RedgeCorrectness::InvalidEdge(i, EdgeCorrectness::AbsentEndpoint(Endpoint::V2));
        }
        if edge.v1_cycle.prev_edge.is_absent() {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CyclePointsToAbsent(Endpoint::V1),
            );
        }
        if edge.v1_cycle.next_edge.is_absent() {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CyclePointsToAbsent(Endpoint::V1),
            );
        }
        if edge.v2_cycle.prev_edge.is_absent() {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CyclePointsToAbsent(Endpoint::V2),
            );
        }
        if edge.v2_cycle.next_edge.is_absent() {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CyclePointsToAbsent(Endpoint::V2),
            );
        }
        if edge.vert_ids[0] == edge.vert_ids[1] {
            return RedgeCorrectness::InvalidEdge(i, EdgeCorrectness::EdgeWithOnlyOnePoint);
        }
        if edge.hedge_id.to_index() >= mesh.hedges_meta.len() && !edge.hedge_id.is_absent() {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::InvalidHedgePointer(edge.hedge_id),
            );
        }

        if edge.hedge_id != HedgeId::ABSENT {
            let hedge_meta = &mesh.hedges_meta[edge.hedge_id.to_index()];
            if hedge_meta.edge_id != edge.id {
                return RedgeCorrectness::InvalidEdge(
                    i,
                    EdgeCorrectness::HedgePointsToDifferentEdge,
                );
            }
        }

        let next_edge_v1 = &mesh.edges_meta[edge.v1_cycle.next_edge.to_index()];
        let prev_edge_v1 = &mesh.edges_meta[edge.v1_cycle.prev_edge.to_index()];
        let next_edge_v2 = &mesh.edges_meta[edge.v2_cycle.next_edge.to_index()];
        let prev_edge_v2 = &mesh.edges_meta[edge.v2_cycle.prev_edge.to_index()];

        let [v1, v2] = edge.vert_ids;
        if next_edge_v1.vert_ids[0] != v1 && next_edge_v1.vert_ids[1] != v1 {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V1, TopologicalRelation::Next, v1),
            );
        }

        if prev_edge_v1.vert_ids[0] != v1 && prev_edge_v1.vert_ids[1] != v1 {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V1, TopologicalRelation::Prior, v1),
            );
        }

        if next_edge_v2.vert_ids[0] != v2 && next_edge_v2.vert_ids[1] != v2 {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V2, TopologicalRelation::Next, v2),
            );
        }

        if prev_edge_v2.vert_ids[0] != v2 && prev_edge_v2.vert_ids[1] != v2 {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V2, TopologicalRelation::Prior, v2),
            );
        }

        if next_edge_v1.cycle(edge.vert_ids[0]).prev_edge != edge.id
            || prev_edge_v1.cycle(edge.vert_ids[0]).next_edge != edge.id
        {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(
                    Endpoint::V1,
                    TopologicalRelation::Next,
                    edge.vert_ids[0],
                ),
            );
        }

        if next_edge_v2.cycle(edge.vert_ids[1]).prev_edge != edge.id
            || prev_edge_v2.cycle(edge.vert_ids[1]).next_edge != edge.id
        {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(
                    Endpoint::V2,
                    TopologicalRelation::Prior,
                    edge.vert_ids[1],
                ),
            );
        }
    }

    // Validate hedges.
    for (i, hedge) in mesh.hedges_meta.iter().enumerate() {
        if !hedge.is_active {
            continue;
        }
        if hedge.id.to_index() != i {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::IdAndIndexMismatch);
        }
        if hedge.face_next_id.is_absent()
            || !mesh.hedges_meta[hedge.face_next_id.to_index()].is_active
        {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::FaceLoopChainIsBroken);
        }
        if hedge.face_prev_id.is_absent()
            || !mesh.hedges_meta[hedge.face_prev_id.to_index()].is_active
        {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::FaceLoopChainIsBroken);
        }
        if hedge.radial_next_id.is_absent()
            || !mesh.hedges_meta[hedge.radial_next_id.to_index()].is_active
        {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::RadialChainIsBroken);
        }
        if hedge.radial_prev_id.is_absent()
            || !mesh.hedges_meta[hedge.radial_prev_id.to_index()].is_active
        {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::RadialChainIsBroken);
        }
        if hedge.source_id.is_absent() || !mesh.verts_meta[hedge.source_id.to_index()].is_active {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::EdgeIsAbsent);
        }
        if hedge.edge_id.is_absent() || !mesh.edges_meta[hedge.edge_id.to_index()].is_active {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::EdgeIsAbsent);
        }

        let mut current = hedge.id;
        let mut iter_count = 0;
        loop {
            if current == HedgeId::ABSENT {
                return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::RadialChainIsBroken);
            }

            let next = mesh.hedges_meta[current.to_index()].radial_next_id;
            if next == HedgeId::ABSENT || !mesh.hedges_meta[next.to_index()].is_active {
                return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::RadialChainIsBroken);
            }

            current = next;
            if current == hedge.id || iter_count > 100 {
                break;
            }
            iter_count += 1;
        }

        if iter_count > 100 {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::LongRadialChain(hedge.id));
        }
    }

    for (i, face) in mesh.faces_meta.iter().enumerate() {
        if !face.is_active {
            continue;
        }
        if face.id.to_index() >= mesh.faces_meta.len() {
            return RedgeCorrectness::InvalidFace(i);
        }
        let hedge_handle = mesh.hedge_handle(face.hedge_id);
        if hedge_handle.face().id() != face.id {
            return RedgeCorrectness::InvalidFace(i);
        }
    }

    for (i, vert) in mesh.verts_meta.iter().enumerate() {
        if !vert.is_active {
            continue;
        }

        let v_handle = mesh.vert_handle(crate::VertId(i));
        for e in v_handle.star_edges() {
            debug_assert!(e.is_active());
        }
    }

    if let Some((v, s1, s2)) = check_edge_vertex_cycles(mesh) {
        return RedgeCorrectness::VertexCyclesDontMatch(v, s1, s2);
    }

    // Verify that iterating through the radial orbit of a edge yields the same number as the total
    // faces that point to that edge.
    let mut edge_face_references = vec![0; mesh.edge_count()];
    // Each face marks all edges it sees (thus increments the expected number of counters by 1).
    for face in &mesh.faces_meta {
        if !face.is_active {
            continue;
        }
        let handle = mesh.face_handle(face.id);
        for h in handle.hedge().face_loop() {
            edge_face_references[h.edge().id().to_index()] += 1;
        }
    }
    // Each radial orbit better see as many faces as the prior loop counted.
    for edge in mesh.meta_edges() {
        if !edge.has_hedge() {
            continue;
        }

        let radial_count = edge.hedge().radial_loop().count();
        assert!(radial_count == edge_face_references[edge.id().to_index()]);
    }

    RedgeCorrectness::Correct
}

/// A regular vertex is a vertex with valence 6.
pub fn count_regular_vertices<R: RedgeContainers>(mesh: &Redge<R>) -> usize {
    mesh.meta_verts()
        .map(|v| (v.star_edges().count() == 6) as usize)
        .sum()
}
/// An isolated vertex is a vertex whose edge pointer is absent.
pub fn count_isolated_vertices<R: RedgeContainers>(mesh: &Redge<R>) -> usize {
    mesh.meta_verts()
        .map(|v| (v.star_edges().count() == 0) as usize)
        .sum()
}

/// Number of inactive vertices O(n).
pub fn count_innactive_vertices<R: RedgeContainers>(mesh: &Redge<R>) -> usize {
    mesh.verts_meta
        .iter()
        .map(|v| (!v.is_active) as usize)
        .sum()
}

/// Compute how many vertices there are for each valence.
pub fn vertex_valence_histogram<R: RedgeContainers>(mesh: &Redge<R>) -> BTreeMap<usize, usize> {
    let mut histogram = BTreeMap::new();
    for v in mesh.meta_verts() {
        let valence = v.star_edges().count();
        let val = histogram.entry(valence).or_default();
        *val += 1;
    }

    histogram
}

/// Topology state of a mesh.
#[derive(Debug, PartialEq, Eq)]
pub enum RedgeManifoldness {
    /// Mesh is two manifold.
    IsManifold,
    /// This is an *error* the redge is broken, the mesh *must* be discarded
    /// and whatever created it debugged.
    IsIncorrect(RedgeCorrectness),
    /// Found a vertex with no incident edges.
    IsolatedVertex(usize),
    /// Found an edge with a numer of incident faces other than 1 or 2.
    NonManifoldEdge(usize, EdgeManifoldness),
}

/// Manifold state of an edge.
#[derive(Debug, PartialEq, Eq)]
pub enum EdgeManifoldness {
    /// Edge is manifold.
    Manifold,
    /// This is an *error* a half edge is pointing to itself as part of the face cycle.
    SelfReferentialFaceCycle,
    /// The `next` field of a half edge in a face cycle points to a half edge whose `prev` is not
    /// itself, or vice versa.
    BrokenFaceLoop,
    /// The `next` field of a half edge in a radial cycle points to a half edge whose `prev` is not
    /// itself, or vice versa.
    BrokenRadialLoop,
}

/// Inspect the mesh and return its manifold state. Some values are hard errors, other
/// merely status reports, see `RedgeManifoldness` for details.
pub fn manifold_state<R: RedgeContainers>(mesh: &Redge<R>) -> RedgeManifoldness {
    match correctness_state(mesh) {
        RedgeCorrectness::Correct => {}
        x => return RedgeManifoldness::IsIncorrect(x),
    }
    for (i, vert) in mesh.verts_meta.iter().enumerate() {
        if !vert.is_active {
            continue;
        }
        if vert.edge_id.is_absent() {
            return RedgeManifoldness::IsolatedVertex(i);
        }
    }

    for (i, hedge) in mesh.hedges_meta.iter().enumerate() {
        if !hedge.is_active {
            continue;
        }

        debug_assert!(hedge.id != HedgeId::ABSENT);

        let hedge_handle = mesh.hedge_handle(hedge.id);

        let radial_neighbours: Vec<_> = hedge_handle.radial_loop().collect();
        if radial_neighbours.len() != 2 && radial_neighbours.len() != 1 {
            return RedgeManifoldness::NonManifoldEdge(i, EdgeManifoldness::BrokenRadialLoop);
        }

        if hedge_handle.face_prev().id() == hedge_handle.id()
            || hedge_handle.face_next().id() == hedge_handle.id()
        {
            return RedgeManifoldness::NonManifoldEdge(
                i,
                EdgeManifoldness::SelfReferentialFaceCycle,
            );
        }

        if hedge_handle.face_prev().id() == hedge_handle.face_next().id() {
            return RedgeManifoldness::NonManifoldEdge(i, EdgeManifoldness::BrokenFaceLoop);
        }
    }

    RedgeManifoldness::IsManifold
}

/// Geometry status of a mesh. Topology correctness is a hard constraint, but geometry correctness
/// is less important.
#[derive(Debug, PartialEq, Eq)]
pub enum GeometryCorrectness {
    /// Geometry is adequate.
    Correct,
    /// Hard error, the topology is broken, the cause *must* be investigated.
    RedgeIsBroken(RedgeCorrectness),
    /// Found two different vertices that share a location in space. May or may not be a problem
    /// depending on the application.
    DuplicatePoints(VertId, VertId),
}

#[derive(Copy, Clone, Debug, Default)]
struct TreePoint<V: InnerSpace<S> + Debug, S: RealField> {
    point: V,
    id: VertId,
    _phantom_data: PhantomData<S>,
}

impl<V: Debug + InnerSpace<S>, S: RealField> PartialEq for TreePoint<V, S> {
    fn eq(&self, other: &Self) -> bool {
        self.point[0] == other.point[0]
            && self.point[1] == other.point[1]
            && self.point[2] == other.point[2]
    }
}

impl<V: InnerSpace<S> + Debug, S: RealField + Signed + Bounded> rstar::Point for TreePoint<V, S> {
    type Scalar = S;
    const DIMENSIONS: usize = 3;

    fn generate(mut generator: impl FnMut(usize) -> Self::Scalar) -> Self {
        let mut v = V::default();
        v[0] = generator(0);
        v[1] = generator(1);
        v[2] = generator(2);

        TreePoint {
            point: v,
            id: VertId::ABSENT,
            _phantom_data: PhantomData,
        }
    }

    fn nth(&self, index: usize) -> Self::Scalar {
        self.point[index]
    }

    fn nth_mut(&mut self, index: usize) -> &mut Self::Scalar {
        &mut self.point[index]
    }
}

/// Inspect the mesh and report whether the geometry is sound. Checks that the topology is also
/// correct.
pub fn validate_geometry_state<R: RedgeContainers, S: RealField + Bounded + Signed>(
    mesh: &Redge<R>,
    epsilon: S,
) -> GeometryCorrectness
where
    VertData<R>: InnerSpace<S>,
{
    match correctness_state(mesh) {
        RedgeCorrectness::Correct => {}
        x => return GeometryCorrectness::RedgeIsBroken(x),
    };

    let mut point_set = RTree::<TreePoint<VertData<R>, S>, _>::new();

    for v in mesh.meta_verts() {
        let point_to_insert = v.data().clone();
        let test_point = TreePoint {
            point: point_to_insert.clone(),
            id: v.id(),
            _phantom_data: PhantomData,
        };
        match point_set.nearest_neighbor(&test_point) {
            Some(tree_point) => {
                if (tree_point.point.clone() - point_to_insert.clone()).norm() <= epsilon {
                    return GeometryCorrectness::DuplicatePoints(tree_point.id, test_point.id);
                } else {
                    point_set.insert(test_point);
                }
            }
            None => {
                point_set.insert(test_point);
            }
        }
    }

    GeometryCorrectness::Correct
}
