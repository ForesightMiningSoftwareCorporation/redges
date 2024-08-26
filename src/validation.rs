// The functions defined here should remain functions, don't refactor them to methods.
// The reason is that, we will likely want to implement trait interfaces
// to abstract meshes. These methods should work for those abstractions.
// If we couple them to the objects it will be harder to do this.

use std::collections::{BTreeMap, BTreeSet};

use crate::{
    container_trait::{PrimitiveContainer, RedgeContainers},
    helpers::check_edge_vertex_cycles,
    EdgeId, Endpoint, HedgeId, Redge, VertId,
};

// TODO: For each function in this file there will be missing checks, add them as
// needed until they are entirely trustworthy.
#[derive(Debug, PartialEq, Eq)]
pub enum RedgeCorrectness {
    Correct,
    MismatchingArrayLengths,
    InvalidVert(usize),
    InvalidEdge(usize, EdgeCorrectness),
    InvalidHedge(usize, HedgeCorrectness),
    InvalidFace(usize),
    VertexCyclesDontMatch(VertId, BTreeSet<EdgeId>, BTreeSet<EdgeId>),
}

#[derive(Debug, PartialEq, Eq)]
pub enum EdgeCorrectness {
    Correct,
    EdgeIsAbsent,
    IdAndIndexMismatch,
    AbsentEndpoint(Endpoint),
    CyclePointsToAbsent(Endpoint),
    CycleIsBroken(Endpoint, VertId),
    EdgeWithOnlyOnePoint,
    InvalidHedgePointer(HedgeId),
    HedgePointsToDifferentEdge,
    RepeatedEndpoint(VertId),
}

#[derive(Debug, PartialEq, Eq)]
pub enum HedgeCorrectness {
    Correct,
    IdAndIndexMismatch,
    FaceLoopChainIsBroken,
    RadialChainIsBroken,
    SourceIsAbsent,
}

/// If this returns `Correct` then it is safe to create handles.
// If the above condition is found to be false ping Makogan and tell him
// to fix asap.
pub fn correctness_state<R: RedgeContainers>(mesh: &Redge<R>) -> RedgeCorrectness {
    if mesh.vert_data.len() != mesh.verts_meta.len() {
        return RedgeCorrectness::MismatchingArrayLengths;
    } else if mesh.edge_data.len() != mesh.edges_meta.len() && mesh.edge_data.len() != 0 {
        return RedgeCorrectness::MismatchingArrayLengths;
    } else if mesh.face_data.len() != mesh.faces_meta.len() && mesh.face_data.len() != 0 {
        return RedgeCorrectness::MismatchingArrayLengths;
    }

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
        if v1 == v2 {
            return RedgeCorrectness::InvalidEdge(i, EdgeCorrectness::RepeatedEndpoint(v1));
        }
        if next_edge_v1.vert_ids[0] != v1 && next_edge_v1.vert_ids[1] != v1 {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V1, v1),
            );
        }

        if prev_edge_v1.vert_ids[0] != v1 && prev_edge_v1.vert_ids[1] != v1 {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V1, v1),
            );
        }

        if next_edge_v2.vert_ids[0] != v2 && next_edge_v2.vert_ids[1] != v2 {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V2, v2),
            );
        }

        if prev_edge_v2.vert_ids[0] != v2 && prev_edge_v2.vert_ids[1] != v2 {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V2, v2),
            );
        }

        if next_edge_v1.cycle(edge.vert_ids[0]).prev_edge != edge.id
            || prev_edge_v1.cycle(edge.vert_ids[0]).next_edge != edge.id
        {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V1, edge.vert_ids[0]),
            );
        }

        if next_edge_v2.cycle(edge.vert_ids[1]).prev_edge != edge.id
            || prev_edge_v2.cycle(edge.vert_ids[1]).next_edge != edge.id
        {
            return RedgeCorrectness::InvalidEdge(
                i,
                EdgeCorrectness::CycleIsBroken(Endpoint::V2, edge.vert_ids[1]),
            );
        }
    }

    for (i, hedge) in mesh.hedges_meta.iter().enumerate() {
        if !hedge.is_active {
            continue;
        }
        if hedge.id.to_index() != i {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::IdAndIndexMismatch);
        }
        if hedge.face_next_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::FaceLoopChainIsBroken);
        }
        if hedge.face_prev_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::FaceLoopChainIsBroken);
        }
        if hedge.radial_next_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::RadialChainIsBroken);
        }
        if hedge.radial_prev_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::RadialChainIsBroken);
        }
        if hedge.source_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i, HedgeCorrectness::SourceIsAbsent);
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

    RedgeCorrectness::Correct
}

/// Returns the number of regular vertices
pub fn count_regular_vertices<R: RedgeContainers>(mesh: &Redge<R>) -> usize {
    let mut regular_verts = 0;
    for v in mesh.meta_verts() {
        let valence = v.star_edges().count();
        regular_verts += (valence == 6) as usize;
    }

    regular_verts
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

#[derive(Debug, PartialEq, Eq)]
pub enum RedgeManifoldness {
    IsManifold,
    IsIncorrect(RedgeCorrectness),
    IsolatedVertex(usize),
    NonManifoldEdge(usize, EdgeManifoldness),
    VertexAssociatedWithMultipleCycles,
}

#[derive(Debug, PartialEq, Eq)]
pub enum EdgeManifoldness {
    Manifold,
    SelfReferentialFaceCycle,
    BrokenFaceLoop,
    BrokenRadialLoop,
}

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

        let radial_neighbours: Vec<_> = hedge_handle.radial_neighbours().collect();
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
