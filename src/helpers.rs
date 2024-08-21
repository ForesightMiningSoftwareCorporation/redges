//! WARNING: These functions can and will break invariants of the Redge
//! use with extreme care.

use std::collections::BTreeSet;

use crate::{
    container_trait::RedgeContainers, EdgeId, Endpoint, FaceId, HedgeId, Redge, StarCycleNode,
    VertId,
};

pub(crate) fn remove_edge_from_cycle<R: RedgeContainers>(
    edge_id: EdgeId,
    endpoint: Endpoint,
    mesh: &mut Redge<R>,
) {
    let active_vertex = match endpoint {
        Endpoint::V1 => mesh.edges_meta[edge_id.to_index()].vert_ids[0],
        Endpoint::V2 => mesh.edges_meta[edge_id.to_index()].vert_ids[1],
    };
    let cycle = mesh.edges_meta[edge_id.to_index()]
        .cycle(active_vertex)
        .clone();

    debug_assert!(cycle.prev_edge != edge_id);
    debug_assert!(cycle.next_edge != edge_id);
    debug_assert!(cycle.next_edge != cycle.prev_edge);

    // Attach the prior and next pointers to each other, thus eliminating
    // all references to the current edge.
    let prior_cycle = mesh.edges_meta[cycle.prev_edge.to_index()].cycle_mut(active_vertex);
    prior_cycle.next_edge = cycle.next_edge;
    debug_assert!(prior_cycle.prev_edge != edge_id);
    debug_assert!(prior_cycle.next_edge != edge_id);
    debug_assert!(prior_cycle.prev_edge != prior_cycle.next_edge);

    let next_cycle = mesh.edges_meta[cycle.next_edge.to_index()].cycle_mut(active_vertex);
    next_cycle.prev_edge = cycle.prev_edge;
    debug_assert!(next_cycle.prev_edge != edge_id);
    debug_assert!(next_cycle.next_edge != edge_id);
    debug_assert!(next_cycle.prev_edge != next_cycle.next_edge);
}

pub(crate) fn join_radial_cycles<R: RedgeContainers>(
    h1: HedgeId,
    h2: HedgeId,
    mesh: &mut Redge<R>,
) {
    let next = mesh.hedges_meta[h1.to_index()].radial_next_id;
    let prev = mesh.hedges_meta[h2.to_index()].radial_prev_id;

    mesh.hedges_meta[next.to_index()].radial_prev_id = prev;
    mesh.hedges_meta[prev.to_index()].radial_next_id = next;

    // Handle self referential cases.
    if next == h1 {
        mesh.hedges_meta[h1.to_index()].radial_prev_id = h2;
        mesh.hedges_meta[h1.to_index()].radial_next_id = h2;
    }
    if prev == h2 {
        mesh.hedges_meta[h2.to_index()].radial_prev_id = h1;
        mesh.hedges_meta[h2.to_index()].radial_next_id = h1;
    }
}

pub(crate) fn remove_hedge_from_radial<R: RedgeContainers>(hedge_id: HedgeId, mesh: &mut Redge<R>) {
    debug_assert!(mesh.hedges_meta[hedge_id.to_index()].is_active);
    let next = mesh.hedges_meta[hedge_id.to_index()].radial_next_id;
    let prev = mesh.hedges_meta[hedge_id.to_index()].radial_prev_id;
    debug_assert!(!next.is_absent());
    debug_assert!(!prev.is_absent());

    mesh.hedges_meta[next.to_index()].radial_prev_id = prev;
    mesh.hedges_meta[prev.to_index()].radial_next_id = next;

    let eid = mesh.hedges_meta[hedge_id.to_index()].edge_id;
    mesh.edges_meta[eid.to_index()].hedge_id = next;
}

pub(crate) fn disable_vert_meta<R: RedgeContainers>(vert_id: VertId, mesh: &mut Redge<R>) {
    let vert = &mut mesh.verts_meta[vert_id.to_index()];
    vert.is_active = false;
    vert.edge_id = EdgeId::new_absent();
}

pub(crate) fn disable_edge_meta<R: RedgeContainers>(edge: EdgeId, mesh: &mut Redge<R>) {
    let edge = &mut mesh.edges_meta[edge.to_index()];

    // Break every single pointer in this edge.
    // Yes its more overhead, but it prevents bugs.
    edge.is_active = false;
    edge.vert_ids = [VertId::new_absent(); 2];
    edge.hedge_id = HedgeId::new_absent();
    edge.v1_cycle = StarCycleNode {
        prev_edge: EdgeId::new_absent(),
        next_edge: EdgeId::new_absent(),
    };
    edge.v2_cycle = StarCycleNode {
        prev_edge: EdgeId::new_absent(),
        next_edge: EdgeId::new_absent(),
    };
}

pub(crate) fn disable_hedge_meta<R: RedgeContainers>(hedge: HedgeId, mesh: &mut Redge<R>) {
    let hedge = &mut mesh.hedges_meta[hedge.to_index()];

    hedge.is_active = false;
    hedge.edge_id = EdgeId::new_absent();
    hedge.face_id = FaceId::new_absent();
    hedge.face_next_id = HedgeId::new_absent();
    hedge.face_prev_id = HedgeId::new_absent();
    hedge.radial_next_id = HedgeId::new_absent();
    hedge.radial_prev_id = HedgeId::new_absent();
    hedge.source_id = VertId::new_absent();
}

pub(crate) fn disable_face_meta<R: RedgeContainers>(face: FaceId, mesh: &mut Redge<R>) {
    let face = &mut mesh.faces_meta[face.to_index()];

    face.is_active = false;
    face.hedge_id = HedgeId::new_absent();
}

fn cycle_endpoint_forward<R: RedgeContainers>(
    mesh: &Redge<R>,
    edge: EdgeId,
    vert: VertId,
) -> Vec<EdgeId> {
    let start = edge;
    let vertex = vert;
    let mut current = start;
    let mut seen = Vec::new();
    loop {
        seen.push(current);
        current = mesh.edges_meta[current.to_index()].cycle(vertex).next_edge;
        if current == start {
            return seen;
        }
    }
}

fn cycle_endpoint_backward<R: RedgeContainers>(
    mesh: &Redge<R>,
    edge: EdgeId,
    vert: VertId,
) -> Vec<EdgeId> {
    let start = edge;
    let vertex = vert;
    let mut current = start;
    let mut seen = Vec::new();
    loop {
        seen.push(current);
        current = mesh.edges_meta[current.to_index()].cycle(vertex).prev_edge;
        if current == start {
            break;
        }
    }

    seen
}

fn set_equality<T>(set1: &BTreeSet<T>, set2: &BTreeSet<T>) -> bool
where
    T: Eq,
{
    let mut equal = true;
    for (a, b) in set1.iter().zip(set2.iter()) {
        equal = equal && (a == b);
    }

    equal && set1.len() == set2.len()
}

pub(crate) fn check_edge_vertex_cycles<R: RedgeContainers>(
    mesh: &Redge<R>,
) -> Option<(VertId, BTreeSet<EdgeId>, BTreeSet<EdgeId>)> {
    let mut vertex_edges_in_cycle = vec![BTreeSet::new(); mesh.verts_meta.len()];
    for edge in mesh.edges_meta.iter() {
        if !edge.is_active {
            continue;
        }

        let mut edge_set_is_fine = |vert_id: VertId| {
            let seen_f = BTreeSet::from_iter(cycle_endpoint_forward(mesh, edge.id, vert_id));
            let seen_b = BTreeSet::from_iter(cycle_endpoint_backward(mesh, edge.id, vert_id));

            if !set_equality(&seen_f, &seen_b) {
                return Some((seen_b, seen_f));
            }

            if vertex_edges_in_cycle[vert_id.to_index()].is_empty() {
                vertex_edges_in_cycle[vert_id.to_index()] = seen_f;
            }

            if !set_equality(&vertex_edges_in_cycle[vert_id.to_index()], &seen_b) {
                return Some((vertex_edges_in_cycle[vert_id.to_index()].clone(), seen_b));
            }

            None
        };

        let [vert1, vert2] = mesh.edges_meta[edge.id.to_index()].vert_ids;

        if let Some((s1, s2)) = edge_set_is_fine(vert1) {
            return Some((vert1, s1, s2));
        }
        if let Some((s1, s2)) = edge_set_is_fine(vert2) {
            return Some((vert1, s1, s2));
        }
    }

    None
}
