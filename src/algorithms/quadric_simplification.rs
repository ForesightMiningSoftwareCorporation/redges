use core::f32;
use core::hash::Hash;
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    ops::Mul,
    usize,
};

use crate::{
    algorithms::{quadrics::Quadric, queue::PQueue},
    container_trait::{
        FaceAttributeGetter, FaceData, PrimitiveContainer, VertData, VertexAttributeGetter,
    },
    edge_handle::EdgeHandle,
    face_handle::{FaceDegeneracies, FaceHandle, FaceMetrics},
    mesh_deleter::MeshDeleter,
    validation::{
        correctness_state, validate_geometry_state, GeometryCorrectness, RedgeCorrectness,
    },
    vert_handle::VertHandle,
    wavefront_loader::ObjData,
    wedge::{self, WedgeDS},
    EdgeId, FaceId, HedgeId, VertId,
};
use linear_isomorphic::prelude::*;
use nalgebra::{constraint, ComplexField, DMatrix, DVector, Vector3, Vector4};
use num::{traits::float::FloatCore, Bounded, Float, Signed};
use num_traits::float::TotalOrder;

use crate::{container_trait::RedgeContainers, Redge};

const EDGE_WEIGHT_PENALTY: f32 = 1_000.0;

#[derive(Debug, PartialEq, Eq)]
pub enum SimplificationStrategy {
    /// Stop simplifying when no valid edges exist for preserving topology (closed meshes remain closed).
    Conservative,
    /// Do whatever it takes to simplify geoemtry for as long as there is geometry to simplify (may introduce boundaries on
    /// otherwise closed meshes).
    Aggressive,
}

#[derive(Debug, PartialEq, Eq)]
pub enum AttributeSimplification {
    NoAttributeSimplification,
    SimplifyAtributes,
}

pub struct QuadricSimplificationConfig {
    pub strategy: SimplificationStrategy,
    pub attribute_simplification: AttributeSimplification,
    pub target_face_count: usize,
}

impl Default for QuadricSimplificationConfig {
    fn default() -> Self {
        Self {
            strategy: SimplificationStrategy::Conservative,
            attribute_simplification: AttributeSimplification::NoAttributeSimplification,
            target_face_count: 10_000,
        }
    }
}

fn tmp_export_to_obj<Vec3: linear_isomorphic::InnerSpace<S>, S: RealField, R: RedgeContainers>(
    vertices: &[Vec3],
    faces: &Vec<FaceData<R>>,
    path: &str,
) -> std::io::Result<()>
where
    FaceData<R>: FaceAttributeGetter<S>,
{
    use std::io::Write;
    let mut file = std::fs::File::create(path)?;
    for v in vertices {
        writeln!(file, "v {} {} {}", v[0], v[1], v[2],)?;
    }

    for face in faces {
        writeln!(file, "vt {} {}", face.attribute(0, 0), face.attribute(0, 1))?;
        writeln!(file, "vt {} {}", face.attribute(1, 0), face.attribute(1, 1))?;
        writeln!(file, "vt {} {}", face.attribute(2, 0), face.attribute(2, 1))?;
    }

    for (i, face) in faces.iter().enumerate() {
        let verts = face.attribute_vertices().to_vec();
        writeln!(
            file,
            "f {}/{} {}/{} {}/{}",
            verts[0].to_index() + 1,
            i * 3 + 1, // vt
            verts[1].to_index() + 1,
            i * 3 + 2, // vt
            verts[2].to_index() + 1,
            i * 3 + 3, // vt
        )?;
    }

    Ok(())
}

// Theory: https://www.cs.cmu.edu/~./garland/Papers/quadrics.pdf
//         https://hhoppe.com/newqem.pdf
/// Uses quadric distances to reduce the total number of edges until a target face count is reached.
pub fn quadric_simplify<S, R, L>(
    mut mesh: Redge<R>,
    config: QuadricSimplificationConfig,
    locked_vertex: L,
) -> (Redge<R>, S)
where
    R: RedgeContainers,
    S: RealField
        + nalgebra::ComplexField
        + FloatCore
        + Mul<VertData<R>, Output = VertData<R>>
        + TotalOrder
        + std::iter::Sum
        + Bounded
        + Signed,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
    L: Fn(VertId, &Redge<R>) -> bool,
{
    println!(
        "radial loop at start {}",
        mesh.edge_handle(EdgeId(267055))
            .hedge()
            .radial_loop()
            .count()
    );

    ObjData::export(
        &vec![
            mesh.edge_handle(EdgeId(267055)).v1().data().clone(),
            mesh.edge_handle(EdgeId(267055)).v2().data().clone(),
        ],
        "edge_at_start.obj",
    );

    // Re-scale geometry for better numerical performance.
    let surface_area: S = mesh.meta_faces().map(|f| f.area()).sum();
    let mean_area = surface_area / S::from(mesh.face_count() as f64).unwrap();

    let scale = S::from(1.0).unwrap() / Float::sqrt(mean_area);
    for i in 0..mesh.vert_data.len() {
        let pos = mesh.vert_data.get(i as u64).clone();
        *mesh.vert_data.get_mut(i as u64) = pos * scale;
    }

    let mut wedges = construct_wedges(&mesh);

    // Start the queue with each edge's cost.
    let mut queue = initialize_queue::<S, R>(&mesh, &wedges);

    let mut deleter = crate::mesh_deleter::MeshDeleter::start_deletion(mesh);

    let mut worst_cost = <S as Float>::min_value();
    let mut dbg = 0;
    while !queue.is_empty() && deleter.active_face_count() > config.target_face_count {
        dbg += 1;
        let (
            cost,
            QueueEdgeData {
                id: eid,
                attributes: optimum,
                wedges_to_merge,
                wedge_order,
            },
        ) = queue.pop().unwrap();

        let edge_handle = deleter.mesh().edge_handle(eid);

        // Skip inactive edges.
        if !edge_handle.is_active() {
            continue;
        }

        // We have selected an edge in an isolated face. So just delete the entire face and move on.
        if edge_handle.hedge().radial_loop().count() == 1 {
            if edge_handle.hedge().face().is_isolated() {
                let fid = edge_handle.hedge().face().id();
                deleter.remove_face(fid);
                continue;
            }

            let total_boundary_edges_in_face = edge_handle
                .hedge()
                .face_loop()
                .fold(0, |acc, h| acc + h.edge().is_boundary() as i32);
            if total_boundary_edges_in_face <= 2 {
                let fid = edge_handle.hedge().face().id();
                deleter.remove_face(fid);
                continue;
            }
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

        let mut geom_optimum = VertData::<R>::default();
        geom_optimum[0] = optimum[0];
        geom_optimum[1] = optimum[1];
        geom_optimum[2] = optimum[2];

        if collapse_would_flip_normal(&edge_handle, &geom_optimum) {
            continue;
        }

        debug_assert!(cost >= S::from(0.).unwrap());
        worst_cost = Float::max(worst_cost, cost);

        let edge_handle = deleter.mesh().edge_handle(eid);
        // Remove all edges that touch the current edge from the queue.
        for e in edge_handle
            .v1()
            .star_edges()
            .chain(edge_handle.v1().link_edges())
            .chain(edge_handle.v2().star_edges())
            .chain(edge_handle.v2().link_edges())
        {
            queue.remove(QueueEdgeData {
                id: e.id(),
                attributes: DVector::default(),
                wedge_order: Vec::new(),
                wedges_to_merge: Vec::new(),
            });
        }

        // Compute the edge collapse and update the new position (and attributes if needed).
        let edge_handle = deleter.mesh.edge_handle(eid);
        wedges.collapse_wedge(&optimum, &edge_handle, &wedges_to_merge, &wedge_order);
        let v1 = edge_handle.v1().id();
        let v2 = edge_handle.v2().id();
        let vid = deleter.collapse_edge(eid);

        let deleted = if vid == v1 {
            v2
        } else if vid == v2 {
            v1
        } else {
            unreachable!("Somehow the new vertex after an edge collapse is not one of the original edge vertices.")
        };

        deleter.update_face_corners(vid);
        sync_wedge_and_redge(vid, deleted, &mut wedges, &mut deleter);

        deleter.mesh().vert_data(vid)[0] = optimum[0];
        deleter.mesh().vert_data(vid)[1] = optimum[1];
        deleter.mesh().vert_data(vid)[2] = optimum[2];

        let faces: Vec<_> = deleter
            .mesh()
            .vert_handle(vid)
            .incident_faces()
            .map(|f| f.id())
            .collect();

        // We just changed the face connectivity by an edge collapse, so
        // we must update the vertex indices of the relevant faces.
        for fid in &faces {
            let data = deleter.mesh.face_data.get_mut(fid.to_index() as u64);
            let vids = data.attribute_vertices_mut().to_vec();

            let face = deleter.mesh.face_handle(*fid);
            let updated: Vec<_> = face.hedge().face_loop().map(|h| h.source().id()).collect();
            // If we find a mismatched vertex index, we need to update the face.
            if let Some(mismatched) = vids.iter().position(|id| !updated.contains(&id)) {
                let data = deleter.mesh.face_data.get_mut(fid.to_index() as u64);
                data.attribute_vertices_mut()[mismatched] = vid;
            }
        }

        let vn = deleter.mesh.vert_handle(vid);
        for e in vn.star_edges().chain(vn.link_edges()) {
            debug_assert!(e.is_active());

            let (mut cost, optimum, wedges_to_merge, wedge_order) =
                edge_cost_with_wedges_final(&e, &wedges);

            let [v1, v2] = e.vertex_ids();
            if locked_vertex(v1, &deleter.mesh) || locked_vertex(v2, &deleter.mesh) {
                cost += S::from(EDGE_WEIGHT_PENALTY).unwrap();
            }

            queue.push(
                cost,
                QueueEdgeData {
                    id: e.id(),
                    attributes: optimum,
                    wedges_to_merge,
                    wedge_order,
                },
            );
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

    let (vert_frag, ..) = deleter.compute_fragmentation_maps();
    let mut res = deleter.end_deletion();
    debug_assert!(correctness_state(&res) == RedgeCorrectness::Correct);

    // Update face connectivity upon termination.
    for face in 0..res.face_count() {
        let f = res.face_data.get_mut(face as u64);
        for v in f.attribute_vertices_mut() {
            *v = VertId(*vert_frag.get(&*v).unwrap());
        }
    }

    (res, worst_cost)
}

fn face_area_after<'r, R: RedgeContainers, S: RealField>(
    face: &FaceHandle<'r, R>,
    vert_id: VertId,
    p: &Vector3<S>,
) -> S
where
    VertData<R>: InnerSpace<S>,
    S: nalgebra::ComplexField,
{
    let mut verts = Vec::new();
    for vert in face.vertices() {
        let id = vert.id();
        let pos = if id == vert_id {
            p.clone()
        } else {
            let d = vert.data().clone();
            Vector3::new(d[0], d[1], d[2])
        };

        verts.push(pos);
    }

    let d1 = verts[1] - verts[0];
    let d2 = verts[2] - verts[0];

    let area = S::from_real(d1.cross(&d2).norm());
    area
}

fn sync_wedge_and_redge<R: RedgeContainers, S: RealField>(
    vid: VertId,
    deleted: VertId,
    wedges: &mut WedgeDS<S>,
    deleter: &mut MeshDeleter<R>,
) where
    FaceData<R>: FaceAttributeGetter<S>,
{
    // Update the wedge with the new topology.
    // All surviving faces that used to have `deleted` as a corner must now have `vid` as that corner instead.
    for face in deleter.mesh.vert_handle(vid).incident_faces() {
        let list = wedges.faces.get_mut(&face.id()).unwrap();
        if list.contains_key(&deleted) {
            let val = list.remove(&deleted).unwrap();
            list.insert(vid, val);
        }
    }

    let face_ids = deleter
        .mesh
        .vert_handle(vid)
        .incident_faces()
        .map(|f| f.id())
        .collect::<Vec<_>>();
    for fid in face_ids {
        let face = deleter.mesh.face_handle(fid);
        let inner_index = face.data().inner_index(vid);
        let wedge = wedges.wedge_from_corner(vid, face.id()).unwrap();
        let fid = face.id();

        for ai in 0..wedge.attributes.len() {
            *deleter.mesh.face_data(fid).attribute_mut(inner_index, ai) = wedge.attributes[ai];
        }
    }
}

struct QueueEdgeData<S> {
    id: EdgeId,
    attributes: DVector<S>,
    /// The wedges that must be merged after an edge collapse.
    wedges_to_merge: Vec<(usize, usize)>,
    /// The order in which the wedge attributes appear in the `attributes` vector.
    /// Needed when reconstructing face attribute values from a wedge.
    wedge_order: Vec<usize>,
}

impl<S> PartialEq for QueueEdgeData<S> {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl<S> Eq for QueueEdgeData<S> {}

impl<S> Hash for QueueEdgeData<S> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

/// For each vertex, compute the quadrics of its incident faces.
fn initialize_queue<S, R>(mesh: &Redge<R>, wedges: &WedgeDS<S>) -> PQueue<QueueEdgeData<S>, S>
where
    R: RedgeContainers,
    S: FloatCore
        + RealField
        + nalgebra::ComplexField
        + Mul<VertData<R>, Output = VertData<R>>
        + TotalOrder,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let mut queue = PQueue::with_capacity(mesh.edge_count());
    for edge in mesh.meta_edges() {
        let (cost, optimum, wedges_to_merge, wedge_order) =
            edge_cost_with_wedges_final(&edge, &wedges);
        queue.push(
            cost,
            QueueEdgeData {
                id: edge.id(),
                attributes: optimum,
                wedges_to_merge,
                wedge_order,
            },
        );
    }

    queue
}

fn edge_cost_with_wedges_final<'r, S, R: RedgeContainers>(
    edge: &EdgeHandle<'r, R>,
    wedges: &WedgeDS<S>,
) -> (S, nalgebra::DVector<S>, Vec<(usize, usize)>, Vec<usize>)
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>> + ComplexField,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let (mut q, mut b, mut final_d, wedges_to_merge, wedge_order) = wedges.wedge_quadric(edge);

    let boundary_penalty = S::from(EDGE_WEIGHT_PENALTY).unwrap();

    // Add a constraint planes for boundary edges.
    if edge.is_boundary() {
        let mut add_phantom_plane = |e: &EdgeHandle<'r, R>| {
            let n = e.hedge().face().unit_normal();
            let e_dir = (e.v2().data().clone() - e.v1().data().clone()).normalized();
            let constraint_normal = e_dir.cross(&n);

            let n = nalgebra::Vector3::new(
                constraint_normal[0],
                constraint_normal[1],
                constraint_normal[2],
            );
            let nn = n * n.transpose();
            let p = e.v1().data().clone();
            let p = nalgebra::Vector3::new(p[0], p[1], p[2]);
            let d = -p.dot(&n);

            let mut top_left = q.view_mut((0, 0), (3, 3));
            let mut top_three = b.rows_mut(0, 3);

            top_left += nn;
            top_three += n * d;
            final_d += d * d;
        };

        // For each boundary edge that touches our active edge, we want
        // the final point to approximate the original boundary conditions as much as possible.
        // So add a phantom plane for each such edge, to nudge the final
        // edge position towards the boundary.
        let boundary_set: BTreeSet<_> = edge
            .v1()
            .star_edges()
            .chain(edge.v2().star_edges())
            .filter(|e| e.is_boundary())
            .map(|e| e.id())
            .collect();
        for eid in boundary_set {
            add_phantom_plane(&edge.redge.edge_handle(eid));
        }
    }

    // When the edge is touching the boundary but it is not in it, we must be careful to
    // not move the boundary.
    let v1_is_fixed = edge.v1().is_in_boundary() && !edge.v2().is_in_boundary();
    let v2_is_fixed = edge.v2().is_in_boundary() && !edge.v1().is_in_boundary();

    let solve = || {
        let det = q.determinant();
        // Safety check to prevent numerical problems when solving the system.
        if Float::abs(det) < S::from(f64::EPSILON * 100.).unwrap() || v1_is_fixed || v2_is_fixed {
            return None;
        }
        if let Some(system) = q.clone().cholesky() {
            Some(system.solve(&-b.clone()))
        } else {
            None
        }
    };

    let (optimal, mut cost) = match solve() {
        Some(optimal) => {
            let cost = (optimal.transpose() * &q * &optimal
                + b.transpose() * &optimal * S::from(2.).unwrap())[0]
                + final_d;

            (optimal, cost)
        }
        // If we could not solve the prior system then we will cheat a little bit and solve an easier one.
        None => {
            let mut target =
                (edge.v1().data().clone() + edge.v2().data().clone()) * S::from(0.5).unwrap();

            if v1_is_fixed {
                target = edge.v1().data().clone();
            }
            if v2_is_fixed {
                target = edge.v2().data().clone();
            }
            let mut res = DVector::zeros(q.nrows());

            let wedges1 = wedges.vertex_wedges(&edge.v1());
            let wedges2 = wedges.vertex_wedges(&edge.v2());
            // Penalize collapsing edges with larger numbers of incident wedges at their endpoints.
            let total_wedges = wedges1.len() + wedges2.len();
            let is_boundary = edge.v1().is_in_boundary() || edge.v2().is_in_boundary();

            res[0] = target[0];
            res[1] = target[1];
            res[2] = target[2];

            let p = nalgebra::Vector3::new(target[0], target[1], target[2]);

            let b_low = b.rows(3, b.nrows() - 3);
            let q_left_low = q.view((3, 0), (q.nrows() - 3, 3));

            let mut diag = q.view((3, 3), (q.nrows() - 3, q.ncols() - 3)).into_owned();
            for i in 0..diag.ncols() {
                diag[(i, i)] = S::from(1.0).unwrap() / diag[(i, i)];
            }

            let s = diag * (-b_low - q_left_low * p);
            assert!(s.len() == res.len() - 3);
            for i in 3..res.len() {
                res[i] = s[i - 3];
            }

            (
                res,
                S::from(
                    (total_wedges as f32 + is_boundary as u8 as f32 * 10.) * EDGE_WEIGHT_PENALTY,
                )
                .unwrap(),
            )
        }
    };

    // See if any of the faces becomes degenerate (i.e. area of 0).
    let p = Vector3::new(optimal[0], optimal[1], optimal[2]);
    let f1 = edge.hedge().face().id();
    let f2 = edge.hedge().radial_next().face().id();
    for face in edge
        .v1()
        .incident_faces()
        .filter(|f| f.id() != f1 && f.id() != f2)
    {
        let post_area = face_area_after(&face, edge.v1().id(), &p);
        if post_area < S::from(f32::EPSILON * 10.).unwrap() {
            cost += S::from(EDGE_WEIGHT_PENALTY * EDGE_WEIGHT_PENALTY).unwrap();
        }
    }
    for face in edge
        .v2()
        .incident_faces()
        .filter(|f| f.id() != f1 && f.id() != f2)
    {
        let post_area = face_area_after(&face, edge.v2().id(), &p);
        if post_area < S::from(f32::EPSILON * 10.).unwrap() {
            cost += S::from(EDGE_WEIGHT_PENALTY * EDGE_WEIGHT_PENALTY).unwrap();
        }
    }

    assert!(cost.is_finite());
    let boundary_test =
        edge.is_boundary() || edge.v1().is_in_boundary() || edge.v2().is_in_boundary();
    (
        cost + S::from(boundary_test as u8 as f32).unwrap() * boundary_penalty,
        optimal,
        wedges_to_merge,
        wedge_order,
    )
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
    let mut adjacent_faces: BTreeSet<_> = edge.v1().incident_faces().map(|f| f.id()).collect();
    adjacent_faces.extend(edge.v2().incident_faces().map(|f| f.id()));

    // Remove from the adjacent faces those that will be deleted by the edge collapse.
    for f in edge.hedge().radial_loop().map(|h| h.face().id()) {
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
            .expect("Internal error: the topology of the redge is broken.")
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

fn construct_wedges<S, R>(mesh: &Redge<R>) -> WedgeDS<S>
where
    R: RedgeContainers,
    S: RealField
        + nalgebra::ComplexField
        + FloatCore
        + Mul<VertData<R>, Output = VertData<R>>
        + TotalOrder
        + std::iter::Sum,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let mut wedges_n = WedgeDS::new();
    for face in mesh.meta_faces() {
        wedges_n.insert_face_attributes(face);
    }

    wedges_n
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

        let (redge, cost) = quadric_simplify(
            redge,
            QuadricSimplificationConfig {
                strategy: SimplificationStrategy::Aggressive,
                attribute_simplification: AttributeSimplification::NoAttributeSimplification,
                target_face_count: 0,
            },
            |_, _| false,
        );

        let (vs, fs, _) = redge.to_face_list();
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
        let (redge, cost) = quadric_simplify(
            redge,
            QuadricSimplificationConfig {
                strategy: SimplificationStrategy::Aggressive,
                attribute_simplification: AttributeSimplification::NoAttributeSimplification,
                target_face_count: 0,
            },
            |_, _| false,
        );
        let duration = start.elapsed();
        println!("Time elapsed in simplification is: {:?}", duration);

        let (vs, fs, _) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_flat_donut.obj");
    }
}
