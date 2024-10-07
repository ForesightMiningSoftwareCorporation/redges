use std::{collections::BTreeSet, ops::Mul};

use crate::{
    algorithms::{quadrics::Quadric, queue::PQueue},
    container_trait::{
        FaceAttributeGetter, FaceData, PrimitiveContainer, VertData, VertexAttributeGetter,
    },
    edge_handle::EdgeHandle,
    face_handle::{FaceDegeneracies, FaceHandle, FaceMetrics},
    mesh_deleter::MeshDeleter,
    validation::{correctness_state, RedgeCorrectness},
    vert_handle::VertHandle,
    wavefront_loader::ObjData,
    EdgeId,
};
use linear_isomorphic::prelude::*;
use nalgebra::{ComplexField, DMatrix, DVector, Vector4};
use num::{traits::float::FloatCore, Float};
use num_traits::float::TotalOrder;

use crate::{container_trait::RedgeContainers, Redge};

const EDGE_WEIGHT_PENALTY: f32 = 100_000.0;

#[derive(Debug, PartialEq, Eq)]
pub enum SimplificationStrategy {
    /// Stop simplifying when no valid edges exist for preserving topology (closed meshes remain closed).
    Conservative,
    /// Do whatever it takes to simplify geoemtry for as long as there is geometry to simplify (may introduce boundaries on
    /// otherwise closed meshes).
    Aggressive,
}

#[derive(Debug, PartialEq, Eq)]
pub enum AttributeBlending {
    NoAttributeBlending,
    BlendAtributes,
}

pub struct QuadricSimplificationConfig {
    pub strategy: SimplificationStrategy,
    pub target_face_count: usize,
}

impl Default for QuadricSimplificationConfig {
    fn default() -> Self {
        Self {
            strategy: SimplificationStrategy::Conservative,
            target_face_count: 10_000,
        }
    }
}

// Theory: https://www.cs.cmu.edu/~./garland/Papers/quadrics.pdf
//         https://hhoppe.com/newqem.pdf
/// Uses quadric distances to reduce the total number of edges until a target face count is reached.
pub fn quadric_simplify<S, R>(mut mesh: Redge<R>, config: QuadricSimplificationConfig) -> Redge<R>
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

    // Start the queue with each edge's cost.
    let mut queue = initialize_queue::<S, R>(&mesh);

    let mut deleter = crate::mesh_deleter::MeshDeleter::start_deletion(mesh);
    while !queue.is_empty() && deleter.active_face_count() > config.target_face_count {
        let (_cost, id) = queue.pop().unwrap();
        let eid = EdgeId(id as usize);

        // Skip inactive edges.
        let edge_handle = deleter.mesh().edge_handle(eid);
        if !edge_handle.is_active() {
            continue;
        }

        // Delete non-manifold edges.
        let edge_handle = deleter.mesh().edge_handle(eid);
        if !edge_handle.has_hedge() {
            deleter.remove_edge(eid);
            continue;
        }

        // Skip bad edges.
        let edge_handle = deleter.mesh().edge_handle(eid);
        if !edge_handle.can_collapse() {
            continue;
        }

        // We will deal witht he boundaries later.
        if edge_handle.is_boundary() {
            continue;
        }

        let (_, optimum) = edge_cost(&edge_handle);
        if collapse_would_flip_normal(&edge_handle, &optimum) {
            continue;
        }

        let edge_handle = deleter.mesh().edge_handle(eid);
        // Remove all edges that touch the current edge from the queue.
        for e in edge_handle
            .v1()
            .star_edges()
            .chain(edge_handle.v1().link_edges())
            .chain(edge_handle.v2().star_edges())
            .chain(edge_handle.v2().link_edges())
        {
            queue.remove(e.id().to_index() as u32);
        }

        let vid = deleter.collapse_edge(eid);
        // Update the position of the collapsed vertex to that wich minimizes the
        // quadric error.
        *deleter.mesh().vert_data(vid) = optimum;

        debug_assert!(correctness_state(&deleter.mesh) == RedgeCorrectness::Correct);

        let vn = deleter.mesh().vert_handle(vid);

        for e in vn.star_edges().chain(vn.link_edges()) {
            debug_assert!(e.is_active());
            let (cost, _new_optimum) = edge_cost(&e);

            queue.push(cost, e.id().to_index() as u32);
        }
    }

    // Doing edge collapse after a certain point is very challenging, as a compromise,
    // if we reach here and we need a smaller mesh, we will just delete faces, if anyone wants
    // to try making an edge collapse that works no matter the situation you have my blessing.
    while deleter.active_face_count() > config.target_face_count
        && config.strategy == SimplificationStrategy::Aggressive
    {
        let face = deleter.mesh.faces_meta.iter().find(|f| f.is_active);
        if let Some(f) = face {
            deleter.remove_face(f.id);
        }
    }

    for i in 0..deleter.mesh.vert_data.len() {
        let pos = deleter.mesh.vert_data.get(i as u64).clone();
        *deleter.mesh.vert_data.get_mut(i as u64) = pos * (S::from(1.0).unwrap() / scale);
    }

    let res = deleter.end_deletion();
    debug_assert!(correctness_state(&res) == RedgeCorrectness::Correct);

    res
}

/// For each vertex, compute the quadrics of its incident faces.
fn initialize_queue<S, R>(mesh: &Redge<R>) -> PQueue<u32, S>
where
    R: RedgeContainers,
    S: FloatCore
        + RealField
        + nalgebra::ComplexField
        + Mul<VertData<R>, Output = VertData<R>>
        + TotalOrder,
    VertData<R>: InnerSpace<S>,
{
    let mut queue = PQueue::with_capacity(mesh.edge_count());
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
    // If this topological operation would break the mesh, add a ludicrous cost.
    if !edge.can_collapse() {
        return (
            S::from(EDGE_WEIGHT_PENALTY * 100_000_000.).unwrap(),
            (edge.v1().data().clone() + edge.v2().data().clone()) * S::from(0.5).unwrap(),
        );
    }
    // Id the edge does not have a hedge, it's an isolated one, we REALLY should be collapsing this one.
    if !edge.has_hedge() {
        return (S::from(S::min_value()).unwrap(), edge.v1().data().clone());
    }
    // When the edge is touching the boundary but it is not in it, we must be careful to
    // not move the boundary.
    if edge.v1().is_in_boundary() && !edge.v2().is_in_boundary() {
        return (
            S::from(EDGE_WEIGHT_PENALTY * EDGE_WEIGHT_PENALTY).unwrap(),
            edge.v1().data().clone(),
        );
    }
    if edge.v2().is_in_boundary() && !edge.v1().is_in_boundary() {
        return (
            S::from(EDGE_WEIGHT_PENALTY * EDGE_WEIGHT_PENALTY).unwrap(),
            edge.v2().data().clone(),
        );
    }

    let vertex_quadric = |vert: &VertHandle<R>| {
        let mut quadric = Quadric::default();

        // For each face incident on the vertex, compute its contribution to the vertex'quadric form.
        // for f in vert.incident_faces() {
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
        let n = edge.hedge().face().unit_normal();
        let e_dir = (edge.v2().data().clone() - edge.v1().data().clone()).normalized();
        let constraint_normal = e_dir.cross(&n);

        let border_quadric = Quadric::from_plane(
            edge.v1().data().clone(),
            constraint_normal.normalized(),
            S::from(constraint_normal.norm() + S::from(EDGE_WEIGHT_PENALTY).unwrap()).unwrap(),
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

fn edge_cost_with_attributes<'r, S, R: RedgeContainers>(
    edge: &EdgeHandle<'r, R>,
) -> (S, nalgebra::DVector<S>)
/*-> (S, VertData<R>)*/
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>> + TotalOrder + ComplexField,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let n = edge.hedge().face().data().attribute_count();
    let mut res_mat = DMatrix::<S>::zeros(n, n);
    let mut res_vec = DVector::<S>::zeros(n);
    let mut res_d = S::from(0.).unwrap();

    for face in edge.v1().incident_faces().chain(edge.v2().incident_faces()) {
        let (m, v, d) = face_quadric_with_attributes(&face);

        res_mat += m;
        res_vec += v;
        res_d += d;
    }

    let optimal = res_mat.clone().lu().solve(&res_vec).unwrap();
    let cost = (optimal.transpose() * res_mat * &optimal
        + res_vec.transpose() * &optimal * S::from(2.).unwrap())[0]
        + res_d * res_d;

    (cost, optimal)
}

// Matrix in page 3, section 4, right column:
// https://hhoppe.com/newqem.pdf
fn face_quadric_with_attributes<'r, S, R: RedgeContainers>(
    face: &FaceHandle<'r, R>,
) -> (DMatrix<S>, DVector<S>, S)
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>> + ComplexField,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let positions: Vec<_> = face
        .vertices()
        .map(|v| {
            [
                v.data().attribute(0),
                v.data().attribute(1),
                v.data().attribute(2),
            ]
        })
        .collect();

    let normal = face.unit_normal();
    let normal = nalgebra::Vector3::new(normal[0], normal[1], normal[2]);

    let p = nalgebra::Vector3::new(positions[0][0], positions[0][1], positions[0][2]);
    let d = normal.dot(&p);
    let dn = normal * d;

    let face_attribute_count = face.data().attribute_count();
    let mat_dim = 3 + face_attribute_count;
    let mut mat = DMatrix::<S>::zeros(mat_dim, mat_dim);
    let mut b = DVector::<S>::zeros(mat_dim);
    b[0] = dn[0];
    b[1] = dn[1];
    b[2] = dn[2];

    let mut corner_mat = normal * normal.transpose();
    let mut cum_d = d * d;
    for i in 0..face_attribute_count {
        let (g, d) = attribute_gradient(&positions, &normal, face, i);
        corner_mat += g * g.transpose();

        mat[(3 + i, 0)] = -g[0];
        mat[(3 + i, 1)] = -g[1];
        mat[(3 + i, 2)] = -g[2];

        mat[(0, 3 + i)] = -g[0];
        mat[(1, 3 + i)] = -g[1];
        mat[(2, 3 + i)] = -g[2];

        mat[(3 + i, 3 + i)] = S::from(1.).unwrap();

        b[0] += g[0] * d;
        b[1] += g[1] * d;
        b[2] += g[2] * d;

        b[3 + i] = -d;

        cum_d += d * d;
    }

    for i in 0..corner_mat.nrows() {
        for j in 0..corner_mat.ncols() {
            mat[(i, j)] += corner_mat[(i, j)];
        }
    }

    (mat, b, cum_d)
}

fn attribute_gradient<'r, S, R: RedgeContainers>(
    positions: &[[S; 3]],
    normal: &nalgebra::Vector3<S>,
    face: &FaceHandle<'r, R>,
    attribute: usize,
) -> (nalgebra::Vector3<S>, S)
where
    FaceData<R>: FaceAttributeGetter<S>,
    VertData<R>: VertexAttributeGetter<S>,
    VertData<R>: InnerSpace<S>,
    S: RealField + ComplexField,
{
    let mut mat = nalgebra::Matrix4::<S>::zeros();
    let mut b = nalgebra::Vector4::<S>::zeros();
    for (i, _vert) in face.vertices().enumerate() {
        debug_assert!(i < 3);

        mat[(i, 0)] = positions[i][0];
        mat[(i, 1)] = positions[i][1];
        mat[(i, 2)] = positions[i][2];
        mat[(i, 3)] = S::from(1.).unwrap();

        b[i] = face.data().attribute(i, attribute);
    }

    mat[(3, 0)] = normal[0];
    mat[(3, 1)] = normal[1];
    mat[(3, 2)] = normal[2];

    let decomp = mat.lu();
    let res = decomp.solve(&b).unwrap();

    (res.fixed_rows::<3>(0).into(), res[3])
}

/// Test if an edge collapse would flip the direction of a face normal.
fn collapse_would_flip_normal<'r, R, S>(edge: &EdgeHandle<'r, R>, new_pos: &VertData<R>) -> bool
where
    R: RedgeContainers,
    S: FloatCore + RealField + nalgebra::ComplexField,
    VertData<R>: InnerSpace<S>,
{
    if !edge.has_hedge() {
        return false;
    }
    // These are all the faces adjacent to the edge.
    let mut adjacent_faces: BTreeSet<_> = edge
        .v1()
        .star_edges()
        .filter(|e| e.has_hedge())
        .flat_map(|e| {
            e.hedge()
                .radial_neighbours()
                .map(|h| h.face().id())
                .collect::<Vec<_>>()
        })
        .collect();
    adjacent_faces.extend(
        edge.v2()
            .star_edges()
            .filter(|e| e.has_hedge())
            .flat_map(|e| {
                e.hedge()
                    .radial_neighbours()
                    .map(|h| h.face().id())
                    .collect::<Vec<_>>()
            }),
    );

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
    fn test_quadric_simplification_closed() {
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

        let redge = quadric_simplify(
            redge,
            QuadricSimplificationConfig {
                strategy: SimplificationStrategy::Aggressive,
                target_face_count: 0,
            },
        );

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_closed.obj");
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
        let redge = quadric_simplify(
            redge,
            QuadricSimplificationConfig {
                strategy: SimplificationStrategy::Aggressive,
                target_face_count: 0,
            },
        );
        let duration = start.elapsed();
        println!("Time elapsed in simplification is: {:?}", duration);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_flat_donut.obj");
    }
}
