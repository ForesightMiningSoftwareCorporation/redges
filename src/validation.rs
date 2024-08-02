// The functions defined here should remain functions, don't refactor them to methods.
// The reason is that, we will likely want to implement trait interfaces
// to abstract meshes. These methods should work for those abstractions.
// If we couple them to the objects it will be harder to do this.

use crate::{edge_handle, hedge_handle, EdgeId, PrimitiveContainer, Redge};

// TODO: For each function in this file there will be missing checks, add them as
// needed until they are entirely trustworthy.
#[derive(Debug, PartialEq, Eq)]
pub enum RedgeCorrectness {
    Correct,
    MismatchingArrayLengths,
    InvalidVert(usize),
    InvalidEdge(usize),
    InvalidHedge(usize),
    InvalidFace(usize),
}

pub fn is_correct<V, E, F>(mesh: &Redge<V, E, F>) -> RedgeCorrectness
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
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
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.vert_ids[0].is_absent() {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.vert_ids[1].is_absent() {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.hedge_id.is_absent() {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.v1_cycle.prev_edge.is_absent() {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.v1_cycle.next_edge.is_absent() {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.v2_cycle.prev_edge.is_absent() {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.v2_cycle.next_edge.is_absent() {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.vert_ids[0] == edge.vert_ids[1] {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if edge.hedge_id.to_index() >= mesh.hedges_meta.len() {
            return RedgeCorrectness::InvalidEdge(i);
        }
        if mesh.hedge_handle(edge.hedge_id).edge().id() != edge.id {
            return RedgeCorrectness::InvalidEdge(i);
        }
    }

    for (i, hedge) in mesh.hedges_meta.iter().enumerate() {
        if !hedge.is_active {
            continue;
        }
        if hedge.id.to_index() != i {
            return RedgeCorrectness::InvalidHedge(i);
        }
        if hedge.face_next_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i);
        }
        if hedge.face_prev_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i);
        }
        if hedge.radial_next_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i);
        }
        if hedge.radial_prev_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i);
        }
        if hedge.source_id.is_absent() {
            return RedgeCorrectness::InvalidHedge(i);
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

    RedgeCorrectness::Correct
}

#[derive(Debug, PartialEq, Eq)]
pub enum RedgeManifoldness {
    IsManifold,
    IsIncorrect(RedgeCorrectness),
    IsolatedVertex(usize),
    NonManifoldEdge(usize),
}

pub fn is_manifold<V, E, F>(mesh: &Redge<V, E, F>) -> RedgeManifoldness
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    match is_correct(mesh) {
        RedgeCorrectness::Correct => {}
        x => return RedgeManifoldness::IsIncorrect(x),
    }

    for (i, vert) in mesh.verts_meta.iter().enumerate() {
        if vert.edge_id.is_absent() {
            return RedgeManifoldness::IsolatedVertex(i);
        }
    }

    for (i, hedge) in mesh.hedges_meta.iter().enumerate() {
        let hedge_handle = mesh.hedge_handle(hedge.id);

        let radial_neighbours: Vec<_> = hedge_handle.radial_neighbours().collect();
        if radial_neighbours.len() != 1 && radial_neighbours.len() != 2 {
            return RedgeManifoldness::NonManifoldEdge(i);
        }
    }

    RedgeManifoldness::IsManifold
}
