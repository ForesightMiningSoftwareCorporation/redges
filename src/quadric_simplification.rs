use std::{
    fmt::Debug,
    fmt::Display,
    ops::{Mul, Neg},
};

use crate::{
    container_trait::{PrimitiveContainer, VertData},
    edge_handle::EdgeHandle,
    face_handle::FaceMetrics,
    p_queue::PQueue,
    validation::{correctness_state, RedgeCorrectness},
    wavefront_loader::ObjData,
    EdgeId,
};
use linear_isomorphic::prelude::*;
use nalgebra::{Matrix4, Vector4};
use num::{traits::float::FloatCore, Float};

use crate::{container_trait::RedgeContainers, Redge};

// Theory: https://www.cs.cmu.edu/~./garland/Papers/quadrics.pdf
/// Uses quadric distances to reduce the total number of edges by a an amount
/// equal to `simplify_count`. For example if `simplify_count` is 1000 then
/// 1000 edges will be collapsed.
pub fn quadric_simplify<S, R>(mesh: Redge<R>, mut simplify_count: usize) -> Redge<R>
where
    R: RedgeContainers,
    S: RealField + nalgebra::ComplexField + FloatCore + Mul<Vector4<S>, Output = Vector4<S>>,
    VertData<R::VertContainer>: InnerSpace<S>,
{
    let (mut quadrics, mut queue) = initialize_vertex_quadrics::<S, _>(&mesh);
    let mut deleter = crate::mesh_deleter::MeshDeleter::start_deletion(mesh);

    while !queue.is_empty() && simplify_count > 0 {
        let (eid, _) = queue.pop().unwrap();

        let handle = deleter.mesh().edge_handle(eid);
        if !handle.is_active() || !handle.can_collapse() {
            continue;
        }

        println!("{}", simplify_count);
        simplify_count -= 1;
        // dbg
        let (vs, fs) = deleter.mesh().to_face_list();
        ObjData::export(&(&vs, &fs), format!("top_{}.obj", simplify_count).as_str());
        //

        let edge_handle = deleter.mesh().edge_handle(eid);

        // Remove all edges that touch the current edge from the queue.
        for e in edge_handle.v1().star_edges() {
            queue.remove(e.id());
        }
        for e in edge_handle.v2().star_edges() {
            queue.remove(e.id());
        }

        let (_cost, point) = compute_edge_weight(&edge_handle, &quadrics);
        let new_quadric =
            quadrics[edge_handle.v1().id().to_index()] + quadrics[edge_handle.v2().id().to_index()];

        // dbg
        let [v1, v2] = edge_handle.vertex_ids();
        let (vs, fs) = deleter.mesh().to_face_list();
        ObjData::export(
            &(&vs, &fs),
            format!("before_remove_{}.obj", simplify_count).as_str(),
        );
        let dbg_edge = vec![
            deleter.mesh().vert_data.get(v1.to_index() as u64).clone(),
            deleter.mesh().vert_data.get(v2.to_index() as u64).clone(),
        ];
        ObjData::export(
            &(&dbg_edge, &vec![[0_usize, 1_usize]]),
            format!("edge_{}.obj", simplify_count).as_str(),
        );
        //

        let v = deleter.collapse_edge(eid);

        // dbg

        let (vs, fs) = deleter.mesh().to_face_list();
        ObjData::export(
            &(&vs, &fs),
            format!("after_remove_{}.obj", simplify_count).as_str(),
        );
        //

        let state = correctness_state(deleter.mesh());
        debug_assert!(state == RedgeCorrectness::Correct, "{:?}", state);
        let vid = deleter.mesh().vert_handle(v).id();
        quadrics[vid.to_index()] = new_quadric;

        // Update the position of the collapsed vertex to that wich minimizes the
        // quadric error.
        let mut pos = VertData::<R::VertContainer>::default();
        debug_assert!(
            Float::is_finite(point.x) && Float::is_finite(point.y) && Float::is_finite(point.z)
        );
        pos[0] = point.x;
        pos[1] = point.y;
        pos[2] = point.z;
        *deleter.mesh().vert_data(vid) = pos;

        // Add all edges touching the collapsed vertex with updated costs.
        for e in deleter.mesh().vert_handle(vid).star_edges() {
            debug_assert!(e.is_active());
            let (cost, _) = compute_edge_weight(&e, &quadrics);
            queue.push(e.id(), cost);
        }
    }

    deleter.end_deletion()
}

/// For each vertex, compute the quadrics around its vertices.
fn initialize_vertex_quadrics<S, R>(mesh: &Redge<R>) -> (Vec<Matrix4<S>>, PQueue<EdgeId, S>)
where
    R: RedgeContainers,
    S: FloatCore + RealField + nalgebra::ComplexField,
    VertData<R::VertContainer>: InnerSpace<S>,
{
    let mut vertex_quadrics = vec![Matrix4::<S>::zeros(); mesh.vert_count()];

    for face in mesh.meta_faces() {
        let normal = face.normal();
        debug_assert!(
            Float::is_finite(normal[0])
                && Float::is_finite(normal[1])
                && Float::is_finite(normal[2])
        );
        let position = face.hedge().source().data().clone();
        let d = position.dot(&normal);

        let q = Vector4::<S>::new(normal[0], normal[1], normal[2], -d);
        let quadric = q * q.transpose();

        for v in face.vertices() {
            vertex_quadrics[v.id().to_index()] += quadric;
        }
    }

    let mut queue = PQueue::new();
    for edge in mesh.meta_edges() {
        let (cost, _) = compute_edge_weight(&edge, &vertex_quadrics);
        queue.push(edge.id(), cost);
    }

    (vertex_quadrics, queue)
}

fn compute_edge_weight<'r, S, R: RedgeContainers>(
    edge: &EdgeHandle<'r, R>,
    quadrics: &Vec<Matrix4<S>>,
) -> (S, Vector4<S>)
where
    S: RealField + nalgebra::ComplexField,
    VertData<R::VertContainer>: VectorSpace<Scalar = S>,
{
    debug_assert!(edge.is_active());
    let cost_matrix = quadrics[edge.v1().id().to_index()] + quadrics[edge.v2().id().to_index()];

    let intermediate = match cost_matrix.try_inverse() {
        None => {
            let half = S::from(0.5).unwrap();
            (edge.v1().data().clone() + edge.v2().data().clone()) * half
        }
        Some(inverse) => {
            //dbg
            let half = S::from(0.5).unwrap();
            (edge.v1().data().clone() + edge.v2().data().clone()) * half

            // let b = Vector4::<S>::new(
            //     S::from(0.).unwrap(),
            //     S::from(0.).unwrap(),
            //     S::from(0.).unwrap(),
            //     S::from(1.).unwrap(),
            // );
            // // Find the point that minimizes the quadric error.
            // let sol = inverse * b;

            // println!("{}", cost_matrix);

            // // Homogenize the coordinates.
            // let mut res = VertData::<R::VertContainer>::default();
            // res[0] = sol.x / sol.w;
            // res[1] = sol.y / sol.w;
            // res[2] = sol.z / sol.w;

            // res
        }
    };

    let new_point = Vector4::<S>::new(
        intermediate[0],
        intermediate[1],
        intermediate[2],
        S::from(1.0).unwrap(),
    );

    (
        -(new_point.transpose() * cost_matrix * new_point)[(0, 0)],
        new_point,
    )
}

#[cfg(test)]
mod tests {
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
        } = ObjData::from_disk_file("assets/loop_cube.obj");

        let vertices: Vec<_> = vertices
            .into_iter()
            .map(|v| Vector3::new(v[0], v[1], v[2]))
            .collect();
        println!("{}", vertices.len());
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

        let redge = quadric_simplify(redge, 30000);

        let (vs, fs) = redge.to_face_list();
        println!("{}", vs.len());
        ObjData::export(&(&vs, &fs), "out/simplified_loop_cube.obj");
    }
}
