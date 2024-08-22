use std::{collections::BTreeSet, hash::Hash, ops::Mul};

use crate::{
    container_trait::{PrimitiveContainer, VertData},
    edge_handle::EdgeHandle,
    face_handle::FaceMetrics,
    p_queue::PQueue,
    quadrics::Quadric,
    wavefront_loader::ObjData,
    EdgeId,
};
use linear_isomorphic::prelude::*;
use num::{traits::float::FloatCore, Float};

use crate::{container_trait::RedgeContainers, Redge};

// Theory: https://www.cs.cmu.edu/~./garland/Papers/quadrics.pdf
/// Uses quadric distances to reduce the total number of edges by a an amount
/// equal to `simplify_count`. For example if `simplify_count` is 1000 then
/// 1000 edges will be collapsed.
pub fn quadric_simplify<S, R>(mut mesh: Redge<R>, mut simplify_count: usize) -> Redge<R>
where
    R: RedgeContainers,
    S: RealField + nalgebra::ComplexField + FloatCore + Mul<VertData<R>, Output = VertData<R>>,
    VertData<R>: InnerSpace<S>,
{
    // let mut surface_area = S::from(0.0).unwrap();
    // for face in mesh.meta_faces() {
    //     surface_area += face.area();
    // }
    // let mean_area = surface_area / S::from(mesh.face_count() as f64).unwrap();
    // for i in 0..mesh.vert_data.len() {
    //     let pos = mesh.vert_data.get(i as u64).clone();
    //     *mesh.vert_data.get_mut(i as u64) = pos * (S::from(1.0).unwrap() / mean_area);
    // }

    let (mut quadrics, mut queue) = initialize_vertex_quadrics::<S, R>(&mesh);

    let mut deleter = crate::mesh_deleter::MeshDeleter::start_deletion(mesh);
    while !queue.is_empty() && simplify_count > 0 {
        let (EdgeOptimum { id: eid, optimum }, _cost) = queue.pop().unwrap();

        let edge_handle = deleter.mesh().edge_handle(eid);
        if !edge_handle.is_active() || !edge_handle.can_collapse() {
            continue;
        }

        if collapse_would_flip_normal(&edge_handle, &optimum) {
            continue;
        }

        simplify_count -= 1;

        let v1 = edge_handle.v1().data().clone();
        let v2 = edge_handle.v2().data().clone();

        let (vs, fs) = deleter.mesh.to_face_list();
        ObjData::export(
            &(&vs, &fs),
            &format!("mesh_{}.obj", simplify_count).to_string(),
        );
        ObjData::export(
            &(&vec![v1.clone(), v2.clone()], &vec![[0, 1]]),
            &format!("edge_{}.obj", simplify_count).to_string(),
        );
        let edge_handle = deleter.mesh().edge_handle(eid);

        // Remove all edges that touch the current edge from the queue.
        for e in edge_handle.v1().star_edges() {
            queue.remove(EdgeOptimum::from_id(e.id()));
        }
        for e in edge_handle.v2().star_edges() {
            queue.remove(EdgeOptimum::from_id(e.id()));
        }

        // Update the quadric.
        let new_quadric = quadrics[edge_handle.v1().id().to_index()].clone()
            + quadrics[edge_handle.v2().id().to_index()].clone();

        let v = deleter.collapse_edge(eid);
        let vid = deleter.mesh().vert_handle(v).id();
        quadrics[vid.to_index()] = new_quadric;

        // Update the position of the collapsed vertex to that wich minimizes the
        // quadric error.
        println!("\n\nbefore {:?}", deleter.mesh().vert_data(vid));
        println!("optimum {:?}", optimum);
        println!("v1 {:?} v2 {:?}", v1, v2);
        *deleter.mesh().vert_data(vid) = optimum;
        println!("after {:?}", deleter.mesh().vert_data(vid));

        // Add all edges touching the collapsed vertex with updated costs.
        for e in deleter.mesh().vert_handle(vid).star_edges() {
            debug_assert!(e.is_active());
            let (cost, new_optimum) = compute_edge_cost(&e, &quadrics);
            queue.push(
                EdgeOptimum {
                    id: e.id(),
                    optimum: new_optimum,
                },
                cost,
            );
        }
    }

    deleter.end_deletion()
}

struct EdgeOptimum<V> {
    id: EdgeId,
    optimum: V,
}

impl<V> PartialEq for EdgeOptimum<V> {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl<V> Eq for EdgeOptimum<V> {}

impl<V> Hash for EdgeOptimum<V> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl<V: Default> EdgeOptimum<V> {
    fn from_id(id: EdgeId) -> Self {
        Self {
            id,
            optimum: V::default(),
        }
    }
}

/// For each vertex, compute the quadrics around its vertices.
fn initialize_vertex_quadrics<S, R>(
    mesh: &Redge<R>,
) -> (
    Vec<Quadric<S, VertData<R>>>,
    PQueue<EdgeOptimum<VertData<R>>, S>,
)
where
    R: RedgeContainers,
    S: FloatCore + RealField + nalgebra::ComplexField + Mul<VertData<R>, Output = VertData<R>>,
    VertData<R>: InnerSpace<S>,
{
    let mut vertex_quadrics = vec![Quadric::default(); mesh.vert_count()];

    for face in mesh.meta_faces() {
        let points: Vec<_> = face.vertices().map(|v| v.data().clone()).take(3).collect();
        for v in face.vertices() {
            vertex_quadrics[v.id().to_index()] +=
                Quadric::from_tri(points[0].clone(), points[1].clone(), points[2].clone());
        }
    }

    let mut queue = PQueue::new();
    for edge in mesh.meta_edges() {
        let (cost, optimum) = compute_edge_cost(&edge, &vertex_quadrics);
        queue.push(
            EdgeOptimum {
                id: edge.id(),
                optimum,
            },
            cost,
        );
    }

    (vertex_quadrics, queue)

    // let mut vertex_quadrics = vec![Matrix4::<S>::zeros(); mesh.vert_count()];

    // for face in mesh.meta_faces() {
    //     let normal = face.normal();
    //     assert!(
    //         Float::is_finite(normal[0])
    //             && Float::is_finite(normal[1])
    //             && Float::is_finite(normal[2])
    //     );
    //     let position = face.hedge().source().data().clone();
    //     let d = -position.dot(&normal);

    //     let q = Vector4::<S>::new(normal[0], normal[1], normal[2], d);
    //     let quadric = q * q.transpose();

    //     for v in face.vertices() {
    //         vertex_quadrics[v.id().to_index()] += quadric;
    //     }
    // }

    // let mut queue = PQueue::new();
    // for edge in mesh.meta_edges() {
    //     let (cost, _) = compute_edge_weight(&edge, &vertex_quadrics);
    //     queue.push(edge.id(), cost);
    // }

    // (vertex_quadrics, queue)
}

fn compute_edge_cost<'r, S, R: RedgeContainers>(
    edge: &EdgeHandle<'r, R>,
    quadrics: &Vec<Quadric<S, VertData<R>>>,
) -> (S, VertData<R>)
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>>,
    VertData<R>: InnerSpace<S>,
{
    let q1 = quadrics[edge.v1().id().to_index()].clone();
    let q2 = quadrics[edge.v2().id().to_index()].clone();

    let q = q1 + q2;

    match q.optimize() {
        Some((s, v)) => (s, v),
        None => {
            let v1 = edge.v1().data().clone();
            let v2 = edge.v2().data().clone();
            let mid = (v1.clone() + v2.clone()) * S::from(0.5).unwrap();

            // let c1 = q.error(v1.clone());
            // let c2 = q.error(v2.clone());
            let cm = q.error(mid.clone());

            // let candidates = [(c1, v1), (c2, v2), (cm, mid)];

            // let best = candidates
            //     .iter()
            //     .min_by(|(a, _), (b, _)| match a.partial_cmp(b) {
            //         Some(x) => x,
            //         None => std::cmp::Ordering::Equal,
            //     })
            //     .unwrap()
            //     .clone();

            // best

            (cm, mid)
        }
    }
}

// fn compute_edge_weight<'r, S, R: RedgeContainers>(
//     edge: &EdgeHandle<'r, R>,
//     quadrics: &Vec<Matrix4<S>>,
// ) -> (S, Vector4<S>)
// where
//     S: RealField + nalgebra::ComplexField,
//     VertData<R>: VectorSpace<Scalar = S>,
// {
//     debug_assert!(edge.is_active());
//     let cost_matrix = quadrics[edge.v1().id().to_index()] + quadrics[edge.v2().id().to_index()];

//     // Alternative point selection when the matrix is singular
//     let pick_safe_point = |edge: &EdgeHandle<R>| {
//         let half = S::from(0.5).unwrap();
//         let mid = (edge.v1().data().clone() + edge.v2().data().clone()) * half;
//         let v1 = edge.v1().data().clone();
//         let v2 = edge.v2().data().clone();

//         let mid = Vector4::new(mid[0], mid[1], mid[2], S::from(1.0).unwrap());
//         let v1 = Vector4::new(v1[0], v1[1], v1[2], S::from(1.0).unwrap());
//         let v2 = Vector4::new(v2[0], v2[1], v2[2], S::from(1.0).unwrap());

//         let c_mid = (mid.transpose() * cost_matrix * mid)[(0, 0)];
//         let c_v1 = (v1.transpose() * cost_matrix * v1)[(0, 0)];
//         let c_v2 = (v2.transpose() * cost_matrix * v2)[(0, 0)];

//         if c_mid < c_v1 && c_mid < c_v2 {
//             return (mid, c_mid);
//         } else if c_v1 < c_v2 {
//             return (v1, c_v1);
//         } else {
//             (v2, c_v2)
//         }
//     };

//     let det = cost_matrix.determinant();
//     let (new_point, cost) = match cost_matrix.try_inverse() {
//         None => pick_safe_point(edge),
//         Some(inverse) => {
//             if Float::abs(det) < S::from(0.0001).unwrap() {
//                 pick_safe_point(edge)
//             } else {
//                 let b = Vector4::<S>::new(
//                     S::from(0.).unwrap(),
//                     S::from(0.).unwrap(),
//                     S::from(0.).unwrap(),
//                     S::from(1.).unwrap(),
//                 );
//                 // Find the point that minimizes the quadric error.
//                 let sol = inverse * b;

//                 // Homogenize the coordinates.
//                 let mut res = VertData::<R>::default();
//                 res[0] = sol.x / sol.w;
//                 res[1] = sol.y / sol.w;
//                 res[2] = sol.z / sol.w;

//                 let new_point = Vector4::<S>::new(res[0], res[1], res[2], S::from(1.0).unwrap());

//                 (
//                     new_point,
//                     (new_point.transpose() * cost_matrix * new_point)[(0, 0)],
//                 )
//             }
//         }
//     };

//     if cost < S::from(0.0).unwrap() {
//         println!("{}", cost);
//     }

//     (-cost, new_point)
// }

/// Test if an edge collapse would flip the direciton of a face normal.
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

        let state = manifold_state(&redge);
        debug_assert!(state == RedgeManifoldness::IsManifold, "{:?}", state);

        let start = Instant::now();
        let redge = quadric_simplify(redge, 250000);
        let duration = start.elapsed();
        println!("Time elapsed in expensive_function() is: {:?}", duration);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_loop_cube.obj");
    }
}
