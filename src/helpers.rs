use crate::{
    container_trait::RedgeContainers, EdgeId, EdgeMetaData, Endpoint, FaceId, HedgeId, Redge,
    StarCycleNode, VertId,
};

/// WARNING: These functions can and will break invariants of the Redge
/// use with extreme care.

pub(crate) fn remove_edge_from_cycle<R: RedgeContainers>(
    edge_id: EdgeId,
    endpoint: Endpoint,
    mesh: &mut Redge<R>,
) {
    let active_vertex = match endpoint {
        Endpoint::V1 => mesh.edges_meta[edge_id.to_index()].vert_ids[0],
        Endpoint::V2 => mesh.edges_meta[edge_id.to_index()].vert_ids[1],
    };
    let cycle = mesh.edges_meta[edge_id.to_index()].cycle(active_vertex);

    // Attach the prior and next poiters to each other, thus elimianting
    // all references to the current edge.
    let prior_cycle = mesh.edges_meta[cycle.prev_edge.to_index()].cycle_mut(active_vertex);
    prior_cycle.next_edge = cycle.next_edge;

    let next_cycle = mesh.edges_meta[cycle.next_edge.to_index()].cycle_mut(active_vertex);
    next_cycle.prev_edge = cycle.prev_edge;
}

pub(crate) fn join_vertex_cycles<R: RedgeContainers>(
    cycle1: StarCycleNode,
    cycle2: StarCycleNode,
    vert_id: VertId,
    mesh: &mut Redge<R>,
) {
    mesh.edges_meta[cycle1.next_edge.to_index()]
        .cycle_mut(vert_id)
        .next_edge = cycle2.prev_edge;
    mesh.edges_meta[cycle2.prev_edge.to_index()]
        .cycle_mut(vert_id)
        .next_edge = cycle1.next_edge;
}

pub(crate) fn join_radial_cycles<R: RedgeContainers>(
    e1: HedgeId,
    e2: HedgeId,
    mesh: &mut Redge<R>,
) {
    let next = mesh.hedges_meta[e1.to_index()].radial_next_id;
    let prev = mesh.hedges_meta[e2.to_index()].radial_prev_id;

    mesh.hedges_meta[next.to_index()].radial_prev_id = prev;
    mesh.hedges_meta[prev.to_index()].radial_next_id = next;
}

pub(crate) fn remove_hedge_from_face<R: RedgeContainers>(hedge_id: HedgeId, mesh: &mut Redge<R>) {
    let next = mesh.hedges_meta[hedge_id.to_index()].face_next_id;
    let prev = mesh.hedges_meta[hedge_id.to_index()].face_prev_id;

    mesh.hedges_meta[next.to_index()].face_prev_id = prev;
    mesh.hedges_meta[prev.to_index()].face_next_id = next;
}

pub(crate) fn remove_hedge_from_radial<R: RedgeContainers>(hedge_id: HedgeId, mesh: &mut Redge<R>) {
    let next = mesh.hedges_meta[hedge_id.to_index()].radial_next_id;
    let prev = mesh.hedges_meta[hedge_id.to_index()].radial_prev_id;

    mesh.hedges_meta[next.to_index()].radial_prev_id = prev;
    mesh.hedges_meta[prev.to_index()].radial_next_id = next;
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
