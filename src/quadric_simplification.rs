use std::{collections::BTreeSet, ops::Mul};

use crate::{
    binary_heap::IndexBinaryHeap,
    container_trait::{PrimitiveContainer, VertData},
    edge_handle::EdgeHandle,
    face_handle::{FaceHandle, FaceMetrics},
    quadrics::Quadric,
    vert_handle::VertHandle,
    EdgeId,
};
use linear_isomorphic::prelude::*;
use num::{traits::float::FloatCore, Float};
use num_traits::float::TotalOrder;

use crate::{container_trait::RedgeContainers, Redge};

const EDGE_WEIGHT: f32 = 100_000.0;

// Theory: https://www.cs.cmu.edu/~./garland/Papers/quadrics.pdf
//         https://hhoppe.com/newqem.pdf
/// Uses quadric distances to reduce the total number of edges by a an amount
/// equal to `simplify_count`. For example if `simplify_count` is 1000 then
/// 1000 edges will be collapsed.
pub fn quadric_simplify<S, R>(mut mesh: Redge<R>, mut simplify_count: usize) -> Redge<R>
where
    R: RedgeContainers,
    S: RealField
        + nalgebra::ComplexField
        + FloatCore
        + Mul<VertData<R>, Output = VertData<R>>
        + TotalOrder
        + std::iter::Sum,
    VertData<R>: InnerSpace<S>,
{
    // Re-scale geometry for better numerical performance.
    let surface_area: S = mesh.meta_faces().map(|f| f.area()).sum();
    let mean_area = surface_area / S::from(mesh.face_count() as f64).unwrap();

    let scale = S::from(1.0).unwrap() / Float::sqrt(mean_area);
    for i in 0..mesh.vert_data.len() {
        let pos = mesh.vert_data.get(i as u64).clone();
        *mesh.vert_data.get_mut(i as u64) = pos * scale;
    }

    let mut queue = initialize_queue::<S, R>(&mesh);

    let mut deleter = crate::mesh_deleter::MeshDeleter::start_deletion(mesh);
    while !queue.is_empty() && simplify_count > 0 {
        let (_cost, id) = queue.pop().unwrap();
        let eid = EdgeId(id as usize);

        let edge_handle = deleter.mesh().edge_handle(eid);
        if !edge_handle.is_active() || !edge_handle.can_collapse() {
            continue;
        }

        let (_, optimum) = edge_cost(&edge_handle);
        if collapse_would_flip_normal(&edge_handle, &optimum) {
            continue;
        }

        simplify_count -= 1;

        let edge_handle = deleter.mesh().edge_handle(eid);

        // Remove all edges that touch the current edge from the queue.
        for e in edge_handle.v1().star_edges() {
            queue.remove(e.id().to_index() as u32);
        }
        for e in edge_handle.v2().star_edges() {
            queue.remove(e.id().to_index() as u32);
        }

        // let v1 = edge_handle.v1().data().clone();
        // let v2 = edge_handle.v2().data().clone();
        // let (vs, fs) = deleter.mesh.to_face_list();
        // ObjData::export(
        //     &(&vs, &fs),
        //     &format!("out/before_mesh_{}.obj", simplify_count).to_string(),
        // );
        // ObjData::export(
        //     &(&vec![v1, v2], &vec![[0, 1]]),
        //     &format!("out/edge_{}.obj", simplify_count).to_string(),
        // );

        let vid = deleter.collapse_edge(eid);

        // Update the position of the collapsed vertex to that wich minimizes the
        // quadric error.
        *deleter.mesh().vert_data(vid) = optimum;

        // let (vs, fs) = deleter.mesh.to_face_list();
        // ObjData::export(
        //     &(&vs, &fs),
        //     &format!("out/after_mesh_{}.obj", simplify_count).to_string(),
        // );

        for e in deleter.mesh().vert_handle(vid).star_edges() {
            debug_assert!(e.is_active());
            let (cost, _new_optimum) = edge_cost(&e);
            queue.push(cost, e.id().to_index() as u32);
        }
    }

    for i in 0..deleter.mesh.vert_data.len() {
        let pos = deleter.mesh.vert_data.get(i as u64).clone();
        *deleter.mesh.vert_data.get_mut(i as u64) = pos * scale;
    }

    deleter.end_deletion()
}

/// For each vertex, compute the quadrics of its incident faces.
fn initialize_queue<S, R>(mesh: &Redge<R>) -> IndexBinaryHeap<S>
where
    R: RedgeContainers,
    S: FloatCore
        + RealField
        + nalgebra::ComplexField
        + Mul<VertData<R>, Output = VertData<R>>
        + TotalOrder,
    VertData<R>: InnerSpace<S>,
{
    let mut queue = IndexBinaryHeap::new();
    for edge in mesh.meta_edges() {
        let (cost, _optimum) = edge_cost(&edge);
        queue.push(cost, edge.id().to_index() as u32);
    }

    queue
}

fn face_quadric<'r, S, R: RedgeContainers>(face: &FaceHandle<'r, R>) -> Quadric<S, VertData<R>>
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>>,
    VertData<R>: InnerSpace<S>,
{
    let points = [
        face.hedge().source().data().clone(),
        face.hedge().face_next().source().data().clone(),
        face.hedge().face_prev().source().data().clone(),
    ];

    Quadric::from_tri(points[0].clone(), points[1].clone(), points[2].clone())
}

fn edge_cost<'r, S, R: RedgeContainers>(edge: &EdgeHandle<'r, R>) -> (S, VertData<R>)
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>> + TotalOrder,
    VertData<R>: InnerSpace<S>,
{
    // When the edge is touching the boundary but it is not in the boundary, we must be careful to
    // not move the boundary.
    if edge.v1().is_in_boundary() && !edge.v2().is_in_boundary() {
        return (
            S::from(EDGE_WEIGHT * EDGE_WEIGHT).unwrap(),
            edge.v1().data().clone(),
        );
    }
    if edge.v2().is_in_boundary() && !edge.v1().is_in_boundary() {
        return (
            S::from(EDGE_WEIGHT * EDGE_WEIGHT).unwrap(),
            edge.v2().data().clone(),
        );
    }

    let vertex_quadric = |vert: &VertHandle<R>| {
        let mut quadric = Quadric::default();

        // For each face incident on the vertex, compute its contribution to the vertex'quadric form.
        for f in vert.incident_faces() {
            quadric += face_quadric(&f);
        }

        quadric
    };

    let q1 = vertex_quadric(&edge.v1());
    let q2 = vertex_quadric(&edge.v2());

    let mut qe = q1 + q2;

    // Add a constraint plane for boundary edges.
    if edge.is_boundary() {
        let n = edge.hedge().face().normal();
        let e_dir = (edge.v2().data().clone() - edge.v1().data().clone()).normalized();
        let constraint_normal = e_dir.cross(&n);

        let border_quadric = Quadric::from_plane(
            edge.v1().data().clone(),
            constraint_normal.normalized(),
            S::from(constraint_normal.norm() * S::from(EDGE_WEIGHT).unwrap()).unwrap(),
        );

        qe += border_quadric;
    }

    let res = match qe.optimize() {
        Some((s, v)) => (s, v),
        None => {
            let v1 = edge.v1().data().clone();
            let v2 = edge.v2().data().clone();
            let mid = (v1.clone() + v2.clone()) * S::from(0.5).unwrap();

            let c1 = qe.error(v1.clone());
            let c2 = qe.error(v2.clone());
            let cm = qe.error(mid.clone());

            let candidates = [(c1, v1), (c2, v2), (cm, mid)];

            let best = candidates
                .iter()
                .min_by(|(a, _), (b, _)| a.total_cmp(b))
                .unwrap()
                .clone();

            best
        }
    };

    res
}

/// Test if an edge collapse would flip the direction of a face normal.
fn collapse_would_flip_normal<'r, R, S>(edge: &EdgeHandle<'r, R>, new_pos: &VertData<R>) -> bool
where
    R: RedgeContainers,
    S: FloatCore + RealField + nalgebra::ComplexField,
    VertData<R>: InnerSpace<S>,
{
    // These are all the faces adjacent to the edge.
    let mut adjacent_faces: BTreeSet<_> = edge
        .v1()
        .star_edges()
        .flat_map(|e| {
            e.hedge()
                .radial_neighbours()
                .map(|h| h.face().id())
                .collect::<Vec<_>>()
        })
        .collect();
    adjacent_faces.extend(edge.v2().star_edges().flat_map(|e| {
        e.hedge()
            .radial_neighbours()
            .map(|h| h.face().id())
            .collect::<Vec<_>>()
    }));

    // Remove from the adjacent faces those that will be deleted by the edge collapse.
    for f in edge.hedge().radial_neighbours().map(|h| h.face().id()) {
        adjacent_faces.remove(&f);
    }

    let get_normal = |p: &VertData<R>, pp: &VertData<R>, pn: &VertData<R>| {
        let e1 = pn.clone() - p.clone();
        let e2 = pp.clone() - p.clone();

        e1.cross(&e2)
    };

    // For each surviving face, see if the new point will flip its orientation.
    let [v1, v2] = edge.vertex_ids();
    for f in adjacent_faces {
        let face = edge.redge.face_handle(f);
        let hedge_at_point = face
            .hedge()
            .face_loop()
            .find(|h| h.source().id() == v1 || h.source().id() == v2)
            .expect("If you are seeing this, the topology of the redge is broken.")
            .id();

        let hedge_at_point = edge.redge.hedge_handle(hedge_at_point);

        let p = hedge_at_point.source().data().clone();
        let pn = hedge_at_point.face_next().source().data().clone();
        let pp = hedge_at_point.face_prev().source().data().clone();

        let current_normal = get_normal(&p, &pp, &pn).normalized();
        let new_normal = get_normal(&new_pos, &pp, &pn).normalized();

        if new_normal.dot(&current_normal) <= S::from(0.0).unwrap() {
            return true;
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use nalgebra::Vector3;

    use crate::validation::{manifold_state, RedgeManifoldness};
    use crate::wavefront_loader::ObjData;

    use super::*;

    #[test]
    fn test_quadric_simplification() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/stanford_dragon.obj");
        // ObjData::from_disk_file("assets/loop_cube.obj");

        let vertices: Vec<_> = vertices
            .into_iter()
            .map(|v| Vector3::new(v[0], v[1], v[2]))
            .collect();

        let redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let start = Instant::now();
        let redge = quadric_simplify(redge, 145000);
        let duration = start.elapsed();
        println!("Time elapsed in expensive_function() is: {:?}", duration);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_loop_cube.obj");
    }

    #[test]
    fn test_quadric_simplification_boundary() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/flat_donut.obj");

        let vertices: Vec<_> = vertices
            .into_iter()
            .map(|v| Vector3::new(v[0], v[1], v[2]))
            .collect();

        let redge = Redge::<(_, _, _)>::new(
            vertices,
            (),
            (),
            vertex_face_indices
                .iter()
                .map(|f| f.iter().map(|&i| i as usize)),
        );

        let start = Instant::now();
        let redge = quadric_simplify(redge, 10000);
        let duration = start.elapsed();
        println!("Time elapsed in expensive_function() is: {:?}", duration);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_flat_donut.obj");
    }
}
