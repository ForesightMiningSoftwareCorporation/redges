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
    EdgeId, FaceId, VertId,
};
use linear_isomorphic::prelude::*;
use nalgebra::{ComplexField, DMatrix, DVector, Vector4};
use num::{traits::float::FloatCore, Bounded, Float, Signed};
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
    // Re-scale geometry for better numerical performance.
    let surface_area: S = mesh.meta_faces().map(|f| f.area()).sum();
    let mean_area = surface_area / S::from(mesh.face_count() as f64).unwrap();

    let scale = S::from(1.0).unwrap() / Float::sqrt(mean_area);
    for i in 0..mesh.vert_data.len() {
        let pos = mesh.vert_data.get(i as u64).clone();
        *mesh.vert_data.get_mut(i as u64) = pos * scale;
    }

    let (vs, _fs, fs_data) = mesh.to_face_list();
    tmp_export_to_obj::<_, _, R>(&vs, &fs_data, "input.obj").unwrap();
    let mut wedges = construct_wedges(&mesh);

    // Start the queue with each edge's cost.
    let mut queue = initialize_queue::<S, R>(&mesh, &config, &wedges);

    let mut deleter = crate::mesh_deleter::MeshDeleter::start_deletion(mesh);
    let (vs, fs, fs_data) = deleter.mesh.to_face_list();
    tmp_export_to_obj::<_, _, R>(&vs, &fs_data, format!("start.obj").as_str()).unwrap();

    let mut dbg = 0;
    let mut worst_cost = <S as Float>::min_value();
    while !queue.is_empty() && deleter.active_face_count() > config.target_face_count {
        dbg += 1;

        let (_cost, id) = queue.pop().unwrap();
        let eid = EdgeId(id as usize);

        let edge_handle = deleter.mesh().edge_handle(eid);
        // Skip inactive edges.
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

        // We will deal with the boundaries later.
        if edge_handle.is_boundary() {
            continue;
        }

        let (cost, optimum) = match config.attribute_simplification {
            AttributeSimplification::NoAttributeSimplification => edge_cost(&edge_handle),
            AttributeSimplification::SimplifyAtributes => {
                let (cost, optimum, _) = edge_cost_with_wedges(&edge_handle, &wedges);

                let mut geom_optimum = VertData::<R>::default();
                geom_optimum[0] = optimum[0];
                geom_optimum[1] = optimum[1];
                geom_optimum[2] = optimum[2];

                (cost, geom_optimum)
            }
        };

        if collapse_would_flip_normal(&edge_handle, &optimum) {
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
            queue.remove(e.id().to_index() as u32);
        }

        let edge_handle = deleter.mesh().edge_handle(eid);

        // Compute the edge collapse and update the new position (and attributes if needed).
        let vid = match config.attribute_simplification {
            AttributeSimplification::NoAttributeSimplification => {
                let (_, optimum) = edge_cost(&edge_handle);
                let vid = deleter.collapse_edge(eid);

                // Update the position of the collapsed vertex to that wich minimizes the
                // quadric error.
                *deleter.mesh().vert_data(vid) = optimum;

                vid
            }
            AttributeSimplification::SimplifyAtributes => {
                let (_, optimum, local_wedge_to_global_wedge) =
                    edge_cost_with_wedges(&edge_handle, &wedges);
                let f1 = edge_handle.hedge().face().id();
                let f2 = edge_handle.hedge().radial_next().face().id();
                let v1 = edge_handle.v1().id();
                let v2 = edge_handle.v2().id();
                // wedges.verify(deleter.mesh());

                let vid = deleter.collapse_edge(eid);
                deleter.update_face_corners(vid);
                wedges.collapse_edge(
                    optimum.clone(),
                    f1,
                    f2,
                    v1,
                    v2,
                    &deleter.mesh.vert_handle(vid),
                    local_wedge_to_global_wedge,
                );

                // debug validation
                for i in 0..deleter.mesh.vert_count() {
                    if !deleter.mesh.face_handle(FaceId(i)).is_active() {
                        continue;
                    }
                    let set1: BTreeSet<_> =
                        fs_data[i].attribute_vertices().iter().copied().collect();
                    let set2: BTreeSet<_> = fs[i].iter().map(|i| VertId(*i as usize)).collect();

                    assert!(set1 == set2);
                }
                tmp_export_to_obj::<_, _, R>(&vs, &fs_data, format!("step_{}.obj", dbg).as_str())
                    .unwrap();
                // std::process::exit(0);
                ObjData::export(&(&vs, &fs), format!("check_{}.obj", dbg).as_str());

                // wedges.verify(deleter.mesh());

                deleter.mesh().vert_data(vid)[0] = optimum[0];
                deleter.mesh().vert_data(vid)[1] = optimum[1];
                deleter.mesh().vert_data(vid)[2] = optimum[2];

                let state = validate_geometry_state(&deleter.mesh, S::from(0.00001).unwrap());
                let (vs, fs, fs_data) = deleter.mesh.to_face_list();
                match state {
                    GeometryCorrectness::DuplicatePoints(v1, v2) => {
                        println!(
                            "duplicates {} {:?} {:?}",
                            dbg,
                            deleter.mesh.vert_handle(v1).data(),
                            deleter.mesh.vert_handle(v2).data(),
                        );
                        let (vs, fs, fs_data) = deleter.mesh.to_face_list();
                        tmp_export_to_obj::<_, _, R>(
                            &vs,
                            &fs_data,
                            format!("very_broken_{}.obj", dbg).as_str(),
                        )
                        .unwrap();
                        ObjData::export(
                            &vec![
                                deleter.mesh.vert_handle(v1).data().clone(),
                                deleter.mesh.vert_handle(v2).data().clone(),
                            ],
                            "points.obj",
                        );
                    }
                    _ => {}
                }
                assert!(state == GeometryCorrectness::Correct, "{:?}", state);

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
                    let updated: Vec<_> =
                        face.hedge().face_loop().map(|h| h.source().id()).collect();
                    // If we find a mismatched vertex index, we need to update the face.
                    if let Some(mismatched) = vids.iter().position(|id| !updated.contains(&id)) {
                        let data = deleter.mesh.face_data.get_mut(fid.to_index() as u64);
                        data.attribute_vertices_mut()[mismatched] = vid;
                    }
                }

                vid
            }
        };

        debug_assert!(correctness_state(&deleter.mesh) == RedgeCorrectness::Correct);

        let vn = deleter.mesh.vert_handle(vid);
        for e in vn.star_edges().chain(vn.link_edges()) {
            debug_assert!(e.is_active());
            let (mut cost, _new_optimum) = edge_cost(&e);
            let [v1, v2] = e.vertex_ids();
            if locked_vertex(v1, &deleter.mesh) || locked_vertex(v2, &deleter.mesh) {
                cost += S::from(EDGE_WEIGHT_PENALTY).unwrap();
            }

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

    // Update the face attributes from the wedge data before returning.
    let active_faces: Vec<_> = deleter.mesh.meta_faces().map(|f| f.id()).collect();
    for fid in active_faces {
        let face_handle = deleter.mesh.face_handle(fid);
        let corners: Vec<_> = face_handle.vertex_ids().collect();
        let wedge_values: Vec<_> = corners
            .iter()
            .map(|vid| wedges.face_corner_wedge(*vid, &face_handle).unwrap())
            .collect();

        for (vindex, val) in wedge_values.iter().enumerate() {
            let face_data = deleter.mesh.face_data(fid);
            for i in 3..val.len() {
                *face_data.attribute_mut(vindex, i - 3) = val[i];
            }
        }
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

/// For each vertex, compute the quadrics of its incident faces.
fn initialize_queue<S, R>(
    mesh: &Redge<R>,
    config: &QuadricSimplificationConfig,
    wedges: &Wedges<S>,
) -> PQueue<u32, S>
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
        match config.attribute_simplification {
            AttributeSimplification::NoAttributeSimplification => {
                let (cost, _optimum) = edge_cost(&edge);
                queue.push(cost, edge.id().to_index() as u32);
            }
            AttributeSimplification::SimplifyAtributes => {
                let (cost, _optimum, _) = edge_cost_with_wedges(&edge, &wedges);
                queue.push(cost, edge.id().to_index() as u32);
            }
        }
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
    // If the edge does not have a hedge, it's an isolated one, we REALLY should be collapsing this one.
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

        // For each face incident on the vertex, compute its contribution to the vertex' quadric form.
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

    qe.optimize().unwrap_or_else(|| {
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
    })
}

fn edge_cost_with_wedges<'r, S, R: RedgeContainers>(
    edge: &EdgeHandle<'r, R>,
    wedges: &Wedges<S>,
) -> (S, nalgebra::DVector<S>, BTreeMap<usize, usize>)
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>> + ComplexField,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let redge = edge.redge;
    let (q, b, d, local_wedge_to_global_wedge) = wedge_quadric(edge.id(), redge, wedges, false);

    let solve = || {
        let det = q.determinant();
        // Safety check to prevent numerical problems when solving the system.
        if Float::abs(det) < S::from(0.1).unwrap() {
            return None;
        }
        q.clone().lu().solve(&-b.clone())
    };

    let (optimal, cost) = match solve() {
        Some(optimal) => {
            let cost = (optimal.transpose() * &q * &optimal
                + b.transpose() * &optimal * S::from(2.).unwrap())[0]
                + d;

            println!("did solve");

            (optimal, cost)
        }
        // If we could not solve the prior system then we will cheat a little bit and solve an easier one.
        None => {
            println!("could not solve");

            let mid = (edge.v1().data().clone() + edge.v2().data().clone()) * S::from(0.5).unwrap();
            let mut res = DVector::zeros(q.nrows());

            res[0] = mid[0];
            res[1] = mid[1];
            res[2] = mid[2];

            let p = nalgebra::Vector3::new(mid[0], mid[1], mid[2]);

            let b_low = b.rows(3, b.nrows() - 3);
            let q_left_low = q.view((3, 0), (q.nrows() - 3, 3));
            let s = -b_low - q_left_low * p;
            assert!(s.len() == res.len() - 3);
            for i in 3..res.len() {
                res[i] = s[i - 3];
            }

            (res, S::from(EDGE_WEIGHT_PENALTY * 100.).unwrap())
        }
    };

    assert!(cost.is_finite());
    (cost, optimal, local_wedge_to_global_wedge)
}

fn wedge_quadric<'r, S, R: RedgeContainers>(
    edge: EdgeId,
    redge: &Redge<R>,
    wedges: &Wedges<S>,
    fix_position: bool,
) -> (DMatrix<S>, DVector<S>, S, BTreeMap<usize, usize>)
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>> + ComplexField,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let edge = redge.edge_handle(edge);
    let v1_id = edge.v1().id();
    let v2_id = edge.v2().id();
    let (wedges_1, faces_1) = wedges.vert_wedges::<R>(edge.v1().id());
    let (_wedges_2, faces_2) = wedges.vert_wedges::<R>(edge.v2().id());

    let woffset = wedges_1.len();
    let mut wedge_faces = BTreeMap::new();

    // The two faces incident on the edge will be removed after the collapse, so ignore them.
    let bad_f1 = edge.hedge().face().id();
    let bad_f2 = edge.hedge().radial_next().face().id();
    for (fid, wid) in faces_1 {
        if fid == bad_f1 || fid == bad_f2 {
            continue;
        }
        let list = wedge_faces.entry(wid).or_insert(Vec::new());
        list.push(redge.face_handle(fid));
    }
    for (fid, wid) in faces_2 {
        if fid == bad_f1 || fid == bad_f2 {
            continue;
        }
        let list = wedge_faces.entry(wid + woffset).or_insert(Vec::new());
        list.push(redge.face_handle(fid));
    }

    // At this point the `wedge_faces` tells us the local association between the kept face sand their wedges.
    let attribute_count = edge.hedge().face().data().attribute_count();
    let n = 3 + wedge_faces.len() * attribute_count;
    let mut final_mat = DMatrix::<S>::zeros(n, n);
    let mut final_b = DVector::<S>::zeros(n);
    let mut cum_d = S::from(0.).unwrap();

    let mut wedge_count = 0;
    let mut local_wedge_to_global_wedge = BTreeMap::<usize, usize>::new();
    for (_, face_list) in wedge_faces {
        let offset = wedge_count * attribute_count;
        assert!(!face_list.is_empty());

        let wid = wedges
            .face_corner_wedge_id(v1_id, &face_list[0])
            .unwrap_or_else(|| wedges.face_corner_wedge_id(v2_id, &face_list[0]).unwrap());
        local_wedge_to_global_wedge.insert(wedge_count, wid);

        for face in face_list {
            if !face.is_active() {
                continue;
            }
            let pos = face.hedge().source().data().clone();
            let p = nalgebra::Vector3::new(pos[0], pos[1], pos[2]);

            // This coefficient will be 0 if the position is fixed, and thus will shut down
            // the geometric coefficients.
            let fixed_pos_coeff = S::from(!fix_position as u8).unwrap();
            let normal = face.unit_normal();
            let mut normal = nalgebra::Vector3::new(normal[0], normal[1], normal[2]);

            // TODO: this is a hack.
            if !normal.x.is_finite() {
                normal.x = S::from(1.0).unwrap();
                normal.y = S::from(0.0).unwrap();
                normal.z = S::from(0.0).unwrap();
            }

            let d = -normal.dot(&p) * fixed_pos_coeff;
            let dn = normal * d * fixed_pos_coeff;

            final_b[0] += dn[0];
            final_b[1] += dn[1];
            final_b[2] += dn[2];

            let mut corner_mat = normal * normal.transpose() * fixed_pos_coeff;
            for i in 0..3 {
                corner_mat[(i, i)] += S::from(fix_position as u8).unwrap();
            }

            cum_d += d * d * fixed_pos_coeff;

            if face.vertices().next().unwrap().id().to_index() == 1276 {
                let (vs, fs, fs_data) = redge.to_face_list();
                tmp_export_to_obj::<_, _, R>(&vs, &fs_data, "very_broken.obj").unwrap();
            }
            let positions: Vec<_> = face
                .vertices()
                .map(|v| {
                    let d = v.data();
                    [d[0], d[1], d[2]]
                })
                .collect();
            let face_attribute_count = face.data().attribute_count();
            let mut face_attribs = [DVector::zeros(0), DVector::zeros(0), DVector::zeros(0)];
            for (i, v) in face.vertex_ids().enumerate() {
                if !face.is_active() {
                    continue;
                }

                face_attribs[i] = wedges.face_corner_wedge(v, &face).unwrap();
            }

            for i in 0..face_attribute_count {
                let (g, d) = attribute_gradient(&positions, &normal, &face_attribs, i);

                corner_mat += g * g.transpose();

                final_mat[(3 + i + offset, 0)] += -g[0];
                final_mat[(3 + i + offset, 1)] += -g[1];
                final_mat[(3 + i + offset, 2)] += -g[2];

                final_mat[(0, 3 + i + offset)] += -g[0];
                final_mat[(1, 3 + i + offset)] += -g[1];
                final_mat[(2, 3 + i + offset)] += -g[2];

                final_mat[(3 + i + offset, 3 + i + offset)] += S::from(1.).unwrap();

                final_b[0] += g[0] * d;
                final_b[1] += g[1] * d;
                final_b[2] += g[2] * d;

                final_b[3 + i + offset] += -d;

                cum_d += d * d;
            }

            for i in 0..corner_mat.nrows() {
                for j in 0..corner_mat.ncols() {
                    final_mat[(i, j)] += corner_mat[(i, j)];
                }
            }
        }

        wedge_count += 1;
    }

    assert!(3 + wedge_count * attribute_count == final_mat.nrows());

    (final_mat, final_b, cum_d, local_wedge_to_global_wedge)
}

fn attribute_gradient<'r, S>(
    positions: &[[S; 3]],
    normal: &nalgebra::Vector3<S>,
    attributes: &[DVector<S>; 3],
    atrribute_id: usize,
) -> (nalgebra::Vector3<S>, S)
where
    S: RealField + ComplexField,
{
    let mut mat = nalgebra::Matrix4::<S>::zeros();
    let mut b = nalgebra::Vector4::<S>::zeros();
    for (i, attribs) in attributes.iter().enumerate() {
        debug_assert!(i < 3);

        mat[(i, 0)] = positions[i][0];
        mat[(i, 1)] = positions[i][1];
        mat[(i, 2)] = positions[i][2];
        mat[(i, 3)] = S::from(1.).unwrap();

        // TODO: Clamp is for uvs, and normals, should not be done in general.
        b[i] = attribs[atrribute_id].clamp(S::from(-1.0).unwrap(), S::from(1.0).unwrap());
    }

    mat[(3, 0)] = normal[0];
    mat[(3, 1)] = normal[1];
    mat[(3, 2)] = normal[2];

    // Trivial solution.
    if b == Vector4::zeros() {
        return (b.fixed_rows::<3>(0).into(), S::from(0.0).unwrap());
    }

    let decomp = mat.lu();
    let res = decomp.solve(&b).expect(format!("{} {}", mat, b).as_str());

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
    let mut adjacent_faces: BTreeSet<_> = edge.v1().incident_faces().map(|f| f.id()).collect();
    adjacent_faces.extend(edge.v2().incident_faces().map(|f| f.id()));

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

struct Wedges<S: RealField> {
    // TODO: this should be extracted from the faces in the redge, so this field should
    // not exist, we are just wasting memory here, but it makes implementation easier for now.
    wedge_values: Vec<DVector<S>>,
    vertices_to_wedges: BTreeMap<VertId, Vec<(usize, FaceId)>>,
}

impl<S: RealField> Wedges<S> {
    fn new() -> Self {
        Self {
            wedge_values: Vec::new(),
            vertices_to_wedges: BTreeMap::new(),
        }
    }

    fn insert_face_attributes<'r, R>(&mut self, vert_id: VertId, face: &FaceHandle<'r, R>)
    where
        R: RedgeContainers,
        VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
        FaceData<R>: FaceAttributeGetter<S>,
        S: ComplexField,
    {
        let wedge_list = self.vertices_to_wedges.entry(vert_id).or_insert(Vec::new());

        let inner_index = face.data().inner_index(vert_id);
        let n = face.data().attribute_count();
        let value = DVector::from_fn(n, |i, _| face.data().attribute(inner_index, i));

        let mut existing_wedge_id = usize::MAX;
        for (wedge_id, _) in wedge_list.iter() {
            let v = &self.wedge_values[*wedge_id];
            if S::from_real((v.clone() - value.clone()).norm()) < S::from(0.1).unwrap() {
                existing_wedge_id = *wedge_id;
                break;
            }
        }

        if existing_wedge_id == usize::MAX {
            self.wedge_values.push(value.clone());
            wedge_list.push((self.wedge_values.len() - 1, face.id()));
        } else {
            wedge_list.push((existing_wedge_id, face.id()));
        }
    }

    fn vert_wedges<'r, R>(&self, vert_id: VertId) -> (Vec<DVector<S>>, Vec<(FaceId, usize)>)
    where
        R: RedgeContainers,
        FaceData<R>: FaceAttributeGetter<S>,
    {
        let map_list = self.vertices_to_wedges.get(&vert_id).unwrap();
        let mut wedge_ids = Vec::new();
        let mut local_associations = Vec::new();

        let mut distinct_count = 0;
        for (wid, fid) in map_list {
            if !wedge_ids.iter().any(|i| *i == *wid) {
                distinct_count += 1;
                wedge_ids.push(*wid);
            }

            local_associations.push((*fid, distinct_count - 1));
        }

        (
            wedge_ids
                .iter()
                .map(|i| self.wedge_values[*i].clone())
                .collect(),
            local_associations,
        )
    }

    // Update wedges so as to get rid of eliminated ones after a collapse.
    fn collapse_edge<'r, R>(
        &mut self,
        optimum: DVector<S>,
        f1: FaceId,
        f2: FaceId,
        v1: VertId,
        v2: VertId,
        vert: &VertHandle<'r, R>,
        local_wedge_to_global_wedge: BTreeMap<usize, usize>,
    ) where
        R: RedgeContainers,
        FaceData<R>: FaceAttributeGetter<S>,
    {
        let deleted = if vert.id() == v1 { v2 } else { v1 };

        // Remove the deleted faces.
        let mut list = self.vertices_to_wedges.get(&v1).unwrap().clone();
        list.retain(|(_, fid)| *fid != f1 && *fid != f2);
        let mut list = self.vertices_to_wedges.get(&v2).unwrap().clone();
        list.retain(|(_, fid)| *fid != f1 && *fid != f2);

        // Delete the removed faces from all wedges.
        let mut list = self.vertices_to_wedges.get(&vert.id()).unwrap().clone();
        list.retain(|(_, fid)| *fid != f1 && *fid != f2);

        for vert in vert.neighbours() {
            let list = self.vertices_to_wedges.get_mut(&vert.id()).unwrap();
            list.retain(|(_, fid)| *fid != f1 && *fid != f2);
        }

        // Get the wedges of the deleted vertex and add them to the new vertex.
        let other_wedges = self.vertices_to_wedges.get(&deleted).unwrap().clone();
        list.extend(
            other_wedges
                .iter()
                .filter(|(_, fid)| *fid != f1 && *fid != f2),
        );

        *self.vertices_to_wedges.get_mut(&vert.id()).unwrap() = list;

        self.vertices_to_wedges.remove(&deleted);

        let attrib_count = vert.edge().hedge().face().data().attribute_count();
        for (local_idx, global_idx) in local_wedge_to_global_wedge {
            let local_attributes =
                optimum.view((3 + local_idx * attrib_count, 0), (attrib_count, 1));
            for i in 0..local_attributes.len() {
                self.wedge_values[global_idx][i] = local_attributes[i];
            }
        }
    }

    fn face_corner_wedge_id<'r, R>(
        &self,
        vert_id: VertId,
        face: &FaceHandle<'r, R>,
    ) -> Option<usize>
    where
        R: RedgeContainers,
        FaceData<R>: FaceAttributeGetter<S>,
    {
        let map_list = self.vertices_to_wedges.get(&vert_id).unwrap();
        map_list
            .iter()
            .find(|(_, fid)| *fid == face.id())
            .map(|(wid, _)| *wid)
    }

    fn face_corner_wedge<'r, R>(
        &self,
        vert_id: VertId,
        face: &FaceHandle<'r, R>,
    ) -> Option<DVector<S>>
    where
        R: RedgeContainers,
        FaceData<R>: FaceAttributeGetter<S>,
    {
        match self.face_corner_wedge_id(vert_id, face) {
            None => None,
            Some(wid) => Some(self.wedge_values[wid].clone()),
        }
    }

    fn wedge_from_id(&self, wid: usize) -> DVector<S> {
        self.wedge_values[wid].clone()
    }

    fn verify<R>(&self, redge: &Redge<R>)
    where
        R: RedgeContainers,
    {
        for vertex in redge.meta_verts() {
            let list = self.vertices_to_wedges.get(&vertex.id()).unwrap().clone();
            for face in vertex.incident_faces() {
                assert!(list.iter().any(|(_, fid)| *fid == face.id()));
            }
        }

        for (vid, list) in &self.vertices_to_wedges {
            for (_, fid) in list {
                assert!(
                    redge.face_handle(*fid).is_active(),
                    "Face {:?}, pointed by {:?}, innactive.",
                    fid,
                    vid
                );
            }
        }
    }
}

fn construct_wedges<S, R>(mesh: &Redge<R>) -> Wedges<S>
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
    let mut wedges = Wedges::new();

    for vert in mesh.meta_verts() {
        for face in vert.incident_faces() {
            wedges.insert_face_attributes(vert.id(), &face);
        }
    }

    wedges
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
                attribute_simplification: AttributeSimplification::NoAttributeSimplification,
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
                attribute_simplification: AttributeSimplification::NoAttributeSimplification,
                target_face_count: 0,
            },
        );
        let duration = start.elapsed();
        println!("Time elapsed in simplification is: {:?}", duration);

        let (vs, fs) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_flat_donut.obj");
    }
}
