//! WARNING: These functions can and will break invariants of the Redge
//! use with extreme care.

use std::{
    collections::{BTreeSet, HashSet},
    ops::Mul,
};

use linear_isomorphic::{InnerSpace, RealField, VectorSpace};

use crate::{
    container_trait::{RedgeContainers, VertData},
    EdgeId, Endpoint, FaceId, HedgeId, Redge, StarCycleNode, VertId,
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

    // This edge should now point to itself since it was removed from the cycle.
    let cycle_mut = mesh.edges_meta[edge_id.to_index()].cycle_mut(active_vertex);
    cycle_mut.next_edge = edge_id;
    cycle_mut.prev_edge = edge_id;

    // Attach the prior and next pointers to each other, thus eliminating
    // all references to the current edge.
    let prior_cycle = mesh.edges_meta[cycle.prev_edge.to_index()].cycle_mut(active_vertex);
    prior_cycle.next_edge = cycle.next_edge;
    let next_cycle = mesh.edges_meta[cycle.next_edge.to_index()].cycle_mut(active_vertex);
    next_cycle.prev_edge = cycle.prev_edge;
}

/// Join two edge endpoint cycles. i.e. Their doubly connected, circular,
/// linked lists are now one. Note, be very conscious of the input order,
/// make sure the `start` parameter is an entirely topologically valid
/// edge.
pub(crate) fn join_vertex_cycles<R: RedgeContainers>(
    head1: EdgeId,
    head2: EdgeId,
    mesh: &mut Redge<R>,
) {
    let pa_meta = &mesh.edges_meta[head1.to_index()];
    let qa_meta = &mesh.edges_meta[head2.to_index()];

    let vertex = pa_meta
        .topology_intersection(qa_meta)
        .expect("Tried to call internal `join_vertex_cycles` with non-overlapping vertices.");

    join_vertex_cycles_at(head1, head2, vertex, mesh);
}

pub(crate) fn join_vertex_cycles_at<R: RedgeContainers>(
    head1: EdgeId,
    head2: EdgeId,
    vertex: VertId,
    mesh: &mut Redge<R>,
) {
    let t1 = mesh.edges_meta[head1.to_index()]
        .cycle(vertex)
        .clone()
        .prev_edge;

    let t2 = mesh.edges_meta[head2.to_index()]
        .cycle(vertex)
        .clone()
        .prev_edge;

    let cycle = mesh.edges_meta[head1.to_index()].cycle_mut(vertex);
    cycle.prev_edge = t2;
    let cycle = mesh.edges_meta[t2.to_index()].cycle_mut(vertex);
    cycle.next_edge = head1;

    let cycle = mesh.edges_meta[t1.to_index()].cycle_mut(vertex);
    cycle.next_edge = head2;
    let cycle = mesh.edges_meta[head2.to_index()].cycle_mut(vertex);
    cycle.prev_edge = t1;
}

pub(crate) fn _collect_forward_cycle<R: RedgeContainers>(
    start: EdgeId,
    vert: VertId,
    mesh: &mut Redge<R>,
) -> Vec<EdgeId> {
    let mut collected = vec![start];
    let mut current = start;
    loop {
        current = mesh.edges_meta[current.to_index()].cycle(vert).next_edge;
        if current == start {
            break;
        }
        collected.push(current);
    }

    collected
}

pub(crate) fn _collect_backward_cycle<R: RedgeContainers>(
    start: EdgeId,
    vert: VertId,
    mesh: &mut Redge<R>,
) -> Vec<EdgeId> {
    let mut collected = vec![start];
    let mut current = start;
    loop {
        current = mesh.edges_meta[current.to_index()].cycle(vert).prev_edge;
        if current == start {
            break;
        }
        collected.push(current);
    }

    collected
}

pub(crate) fn join_radial_cycles<R: RedgeContainers>(
    head1: HedgeId,
    head2: HedgeId,
    mesh: &mut Redge<R>,
) {
    let t1 = mesh.hedges_meta[head1.to_index()].radial_prev_id;
    let t2 = mesh.hedges_meta[head2.to_index()].radial_prev_id;

    mesh.hedges_meta[head1.to_index()].radial_prev_id = t2;
    mesh.hedges_meta[t2.to_index()].radial_next_id = head1;

    mesh.hedges_meta[t1.to_index()].radial_next_id = head2;
    mesh.hedges_meta[head2.to_index()].radial_prev_id = t1;
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
    // If the next and prev point to this hedge, this hedge must have been the last one orbiting the edge.
    mesh.edges_meta[eid.to_index()].hedge_id = if next == hedge_id {
        debug_assert!(prev == hedge_id);
        HedgeId::ABSENT
    } else {
        next
    };
}

pub(crate) fn disable_vert_meta<R: RedgeContainers>(vert_id: VertId, mesh: &mut Redge<R>) {
    if vert_id == VertId::ABSENT {
        return;
    }

    let vert = &mut mesh.verts_meta[vert_id.to_index()];
    vert.is_active = false;
    vert.edge_id = EdgeId::ABSENT;
}

pub(crate) fn disable_edge_meta<R: RedgeContainers>(edge: EdgeId, mesh: &mut Redge<R>) {
    if edge == EdgeId::ABSENT {
        return;
    }

    let edge = &mut mesh.edges_meta[edge.to_index()];

    // Break every single pointer in this edge.
    // Yes its more overhead, but it prevents bugs.
    edge.is_active = false;
    edge.vert_ids = [VertId::ABSENT; 2];
    edge.hedge_id = HedgeId::ABSENT;
    edge.v1_cycle = StarCycleNode {
        prev_edge: EdgeId::ABSENT,
        next_edge: EdgeId::ABSENT,
    };
    edge.v2_cycle = StarCycleNode {
        prev_edge: EdgeId::ABSENT,
        next_edge: EdgeId::ABSENT,
    };
}

pub(crate) fn disable_hedge_meta<R: RedgeContainers>(hedge: HedgeId, mesh: &mut Redge<R>) {
    if hedge == HedgeId::ABSENT {
        return;
    }
    let hedge = &mut mesh.hedges_meta[hedge.to_index()];

    hedge.is_active = false;
    hedge.edge_id = EdgeId::ABSENT;
    hedge.face_id = FaceId::ABSENT;
    hedge.face_next_id = HedgeId::ABSENT;
    hedge.face_prev_id = HedgeId::ABSENT;
    hedge.radial_next_id = HedgeId::ABSENT;
    hedge.radial_prev_id = HedgeId::ABSENT;
    hedge.source_id = VertId::ABSENT;
}

pub(crate) fn disable_face_meta<R: RedgeContainers>(face: FaceId, mesh: &mut Redge<R>) {
    if face == FaceId::ABSENT {
        return;
    }

    let face = &mut mesh.faces_meta[face.to_index()];

    face.is_active = false;
    face.hedge_id = HedgeId::ABSENT;
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
            return seen;
        }
    }
}

fn set_equality<T>(set1: &BTreeSet<T>, set2: &BTreeSet<T>) -> bool
where
    T: Eq,
{
    set1 == set2
}

pub(crate) fn check_edge_vertex_cycles<R: RedgeContainers>(
    mesh: &Redge<R>,
) -> Option<(VertId, BTreeSet<EdgeId>, BTreeSet<EdgeId>)> {
    let mut vertex_edges_in_cycle = vec![BTreeSet::new(); mesh.verts_meta.len()];
    for edge in mesh.edges_meta.iter() {
        if !edge.is_active {
            continue;
        }

        let mut check_edge_set = |vert_id: VertId| {
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

        if let Some((s1, s2)) = check_edge_set(vert1) {
            return Some((vert1, s1, s2));
        }
        if let Some((s1, s2)) = check_edge_set(vert2) {
            return Some((vert1, s1, s2));
        }
    }

    None
}

pub(crate) fn link_hedges_in_face<R: RedgeContainers>(
    h1: HedgeId,
    h2: HedgeId,
    mesh: &mut Redge<R>,
) {
    mesh.hedges_meta[h1.to_index()].face_next_id = h2;
    mesh.hedges_meta[h2.to_index()].face_prev_id = h1;
}

pub(crate) fn link_face<R: RedgeContainers>(hedges: &[HedgeId], face: FaceId, mesh: &mut Redge<R>) {
    for i in 0..hedges.len() {
        let h1 = hedges[i];
        let h2 = hedges[(i + 1) % hedges.len()];
        link_hedges_in_face(h1, h2, mesh);

        mesh.hedges_meta[h1.to_index()].face_id = face;
    }

    mesh.faces_meta[face.to_index()].hedge_id = hedges[0];
}

pub(crate) fn hedge_collapse<R: RedgeContainers>(hedge_id: HedgeId, mesh: &mut Redge<R>) {
    let h0 = hedge_id;
    let hn = mesh.hedges_meta[h0.to_index()].face_next_id;
    let hp = mesh.hedges_meta[h0.to_index()].face_prev_id;

    let f = mesh.hedges_meta[h0.to_index()].face_id;

    remove_hedge_from_radial(h0, mesh);
    disable_hedge_meta(h0, mesh);

    mesh.hedges_meta[hn.to_index()].face_prev_id = hp;
    mesh.hedges_meta[hp.to_index()].face_next_id = hn;

    mesh.faces_meta[f.to_index()].hedge_id = hn;
}

pub(crate) fn digon_to_edge<R: RedgeContainers>(face_id: FaceId, mesh: &mut Redge<R>) {
    let handle = mesh.face_handle(face_id);
    debug_assert!(handle.hedge().face_loop().count() == 2);

    let h1 = handle.hedge().id();
    let h2 = handle.hedge().face_next().id();

    let e1 = handle.hedge().edge().id();
    let e2 = handle.hedge().face_next().edge().id();

    let v1 = handle.hedge().source().id();
    let v2 = handle.hedge().face_next().source().id();

    let h1_safe = handle
        .hedge()
        .radial_loop()
        .find(|h| h.id() != h1)
        .map(|h| h.id())
        .unwrap_or(HedgeId::ABSENT);
    let h2_safe = handle
        .hedge()
        .face_next()
        .radial_loop()
        .find(|h| h.id() != h2)
        .map(|h| h.id())
        .unwrap_or(HedgeId::ABSENT);

    let h2_radials: Vec<_> = mesh
        .hedge_handle(h2)
        .radial_loop()
        .map(|h| h.id())
        .filter(|h| *h != h2)
        .collect();

    // Joining first matters, maintain this order of operations.
    // (It's simpler to remove from a large linked list than to add to a small one).
    if h1_safe != HedgeId::ABSENT && h2_safe != HedgeId::ABSENT {
        join_radial_cycles(h1_safe, h2_safe, mesh);
    }

    if h1_safe != HedgeId::ABSENT {
        remove_hedge_from_radial(h1, mesh);
    }
    if h2_safe != HedgeId::ABSENT {
        remove_hedge_from_radial(h2, mesh);
    }

    for hid in h2_radials {
        mesh.hedges_meta[hid.to_index()].edge_id = e1;
    }

    // Make sure the vertices remain valid.
    mesh.verts_meta[v1.to_index()].edge_id = e1;
    mesh.verts_meta[v2.to_index()].edge_id = e1;

    mesh.edges_meta[e1.to_index()].hedge_id = if h1_safe != HedgeId::ABSENT {
        h1_safe
    } else {
        h2_safe
    };

    remove_edge_from_cycle(e2, Endpoint::V1, mesh);
    remove_edge_from_cycle(e2, Endpoint::V2, mesh);

    disable_edge_meta(e2, mesh);
    disable_hedge_meta(h1, mesh);
    disable_hedge_meta(h2, mesh);
    disable_face_meta(face_id, mesh);
}

pub(crate) fn digon_holes_to_edge<R: RedgeContainers>(mut edges: Vec<EdgeId>, mesh: &mut Redge<R>) {
    let canonical_edge = edges.pop().unwrap();

    let handle = mesh.edge_handle(canonical_edge);
    let v1 = handle.v1().id();
    let v2 = handle.v2().id();

    for other_id in edges {
        let handle = mesh.edge_handle(canonical_edge);
        let other = mesh.edge_handle(other_id);

        let h1_safe = if handle.has_hedge() {
            handle.hedge().id()
        } else {
            HedgeId::ABSENT
        };
        let h2_safe = if other.has_hedge() {
            other.hedge().id()
        } else {
            HedgeId::ABSENT
        };
        let other = other.id();

        // Joining first matters, maintain this order of operations.
        // (It's simpler to remove from a large linked list than to add to a small one).
        if h1_safe != HedgeId::ABSENT && h2_safe != HedgeId::ABSENT {
            join_radial_cycles(h1_safe, h2_safe, mesh);
        }

        let h2_radials: Vec<_> = if h2_safe != HedgeId::ABSENT {
            mesh.hedge_handle(h2_safe)
                .radial_loop()
                .map(|h| h.id())
                .collect()
        } else {
            Vec::new()
        };

        for hid in h2_radials {
            mesh.hedges_meta[hid.to_index()].edge_id = canonical_edge;
        }

        // Make sure the vertices remain valid.
        mesh.verts_meta[v1.to_index()].edge_id = canonical_edge;
        mesh.verts_meta[v2.to_index()].edge_id = canonical_edge;

        mesh.edges_meta[canonical_edge.to_index()].hedge_id = if h1_safe != HedgeId::ABSENT {
            h1_safe
        } else {
            h2_safe
        };

        remove_edge_from_cycle(other, Endpoint::V1, mesh);
        remove_edge_from_cycle(other, Endpoint::V2, mesh);
        disable_edge_meta(other, mesh);
    }
}

/// Returns, in this order, the two half edges along the direction of the input hedge and the new transversal edge.
pub(crate) fn split_hedge<R: RedgeContainers>(
    new_vert_id: VertId,
    hedge_id: HedgeId,
    mesh: &mut Redge<R>,
) -> ([HedgeId; 2], EdgeId) {
    let handle = mesh.hedge_handle(hedge_id);

    let h1 = hedge_id;
    let h2 = handle.face_next().id();
    let h3 = handle.face_prev().id();

    let v1 = handle.face_next().source().id();
    let v2 = handle.source().id();
    let t = handle.face_prev().source().id();

    let edge_at_t = handle.face_prev().source().edge().id();

    let f = handle.face().id();
    let f_data = handle.face().data().clone();

    let vn = new_vert_id;
    let data = handle.edge().data().clone();

    let n1 = mesh.add_hedge(t);
    let n2 = mesh.add_hedge(vn);
    let n3 = mesh.add_hedge(v2);

    let f_n = mesh.add_face(f_data);

    // Create the new face loop.
    link_hedges_in_face(n3, n2, mesh);
    link_hedges_in_face(n2, h3, mesh);
    link_hedges_in_face(h3, n3, mesh);

    // Update the source of the new face hedges.
    mesh.hedges_meta[n3.to_index()].source_id = v2;
    mesh.hedges_meta[n2.to_index()].source_id = vn;
    mesh.hedges_meta[h3.to_index()].source_id = t;

    // New face hedges must point to the new face.
    mesh.hedges_meta[n3.to_index()].face_id = f_n;
    mesh.hedges_meta[n2.to_index()].face_id = f_n;
    mesh.hedges_meta[h3.to_index()].face_id = f_n;

    // Create the old face loop.
    link_hedges_in_face(h1, h2, mesh);
    link_hedges_in_face(h2, n1, mesh);
    link_hedges_in_face(n1, h1, mesh);

    // Update the source of the old face hedges.
    mesh.hedges_meta[h1.to_index()].source_id = vn;
    mesh.hedges_meta[h2.to_index()].source_id = v1;
    mesh.hedges_meta[n1.to_index()].source_id = t;

    // Old face hedges must point to the old face.
    mesh.hedges_meta[h1.to_index()].face_id = f;
    mesh.hedges_meta[h2.to_index()].face_id = f;
    mesh.hedges_meta[n1.to_index()].face_id = f;

    // update the faces to have a correct hedge.
    mesh.faces_meta[f_n.to_index()].hedge_id = n2;
    mesh.faces_meta[f.to_index()].hedge_id = n1;

    join_radial_cycles(n2, n1, mesh);

    // Transversal edge.
    let et = mesh.add_edge([vn, t], data);
    mesh.edges_meta[et.to_index()].hedge_id = n1;
    mesh.hedges_meta[n1.to_index()].edge_id = et;
    mesh.hedges_meta[n2.to_index()].edge_id = et;

    mesh.edges_meta[et.to_index()].hedge_id = n1;

    join_vertex_cycles(edge_at_t, et, mesh);

    ([n3, h1], et)
}
