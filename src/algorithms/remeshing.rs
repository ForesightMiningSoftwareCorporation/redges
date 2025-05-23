//! Remeshing algorithms, used to increase resolution of a mesh while simultaneously trying to improve
//! discrete topology.
use core::f64::consts::PI;
use core::fmt::Debug;
use linear_isomorphic::{InnerSpace, RealField, VectorSpace};
use ordered_float::FloatCore;
use std::{collections::HashSet, ops::Mul};

use super::pqueue::PQueue;
use crate::container_trait::PrimitiveContainer;
use crate::{
    container_trait::{RedgeContainers, VertData},
    edge_handle::EdgeHandle,
    validation::{correctness_state, RedgeCorrectness},
    vert_handle::VertHandle,
    EdgeId, FaceId, Redge, VertId,
};
/// Input argument for remeshing.
pub struct RemeshingContext<S: RealField + FloatCore, R: RedgeContainers> {
    /// Priority queue.
    pub queue: PQueue<EdgeId, S>,
    /// Mesh that will be remeshed.
    pub mesh: Redge<R>,
    /// Vertices that have been modified since the last call to remeshing.
    pub modified_vertices: Vec<VertId>,
    /// Faces that have been modified since the last call to remeshing.
    pub modified_faces: Vec<FaceId>,
}

/// Input argument for remeshing without edge collapsing.
pub struct RemeshingParametersWithoutCollapse<S>
where
    S: Debug + Clone + Default + num::Float,
{
    /// Don't allow triangles with angles smaller than this value.
    pub minimum_triangle_angle: S,
    /// Edges whose incident triangles
    pub dihedral_angle_tolerance: S,
    /// How many new vertices should we try to add before returning.
    pub target_additional_vertices: usize,
}

impl<S> Default for RemeshingParametersWithoutCollapse<S>
where
    S: Debug + Clone + Default + num::Float,
{
    fn default() -> Self {
        Self {
            minimum_triangle_angle: S::from(2.0 * PI / 9.0).unwrap(),
            dihedral_angle_tolerance: S::from(2.0 * PI / 4.0).unwrap(),
            target_additional_vertices: 4000,
        }
    }
}

impl<S, R> RemeshingContext<S, R>
where
    R: RedgeContainers,
    S: RealField + num::traits::float::FloatCore,
    VertData<R>: InnerSpace<S>,
{
    /// Construct context from mesh.
    pub fn from_mesh(mesh: Redge<R>) -> Self {
        debug_assert!(correctness_state(&mesh) == RedgeCorrectness::Correct);
        let pq = PQueue::from_iterator(mesh.meta_edges().map(|e| {
            (
                e.id(),
                (e.v1().data().clone() - e.v2().data().clone()).norm(),
            )
        }));

        Self {
            queue: pq,
            mesh,
            modified_vertices: Vec::new(),
            modified_faces: Vec::new(),
        }
    }
}

/// Strictly increase the number of vertices in a mesh via subdivision. Note that it is possible
/// that if this method returns 0, that calling the method again will produce more edges due
/// to equalization and relaxation.
/// Currently, reprojection doesn't do anything.
///
/// Returns the number of added elements (could be less than `parameters.target_additional_vertices`).
pub fn incremental_refinement_with_context<R: RedgeContainers, S, L, P>(
    context: &mut RemeshingContext<S, R>,
    parameters: RemeshingParametersWithoutCollapse<S>,
    adaptive_target_length: L,
    _reproject: P,
) -> usize
where
    S: RealField + num::traits::float::FloatCore + Mul<VertData<R>, Output = VertData<R>>,
    VertData<R>: InnerSpace<S>,
    L: Fn(usize, usize) -> S,
    P: Fn(&VertData<R>) -> VertData<R>,
{
    debug_assert!(correctness_state(&context.mesh) == RedgeCorrectness::Correct);

    let (mut feat_verts, mut feat_edges, _corners) =
        detect_features::<S, R>(&mut context.mesh, parameters.dihedral_angle_tolerance);

    let high = S::from(4.0 / 3.0).unwrap();

    let input_verts = context.mesh.vert_count();
    while context.mesh.vert_count() < input_verts + parameters.target_additional_vertices {
        if context.queue.is_empty() {
            let pq = PQueue::from_iterator(context.mesh.meta_edges().map(|e| {
                (
                    e.id(),
                    (e.v1().data().clone() - e.v2().data().clone()).norm(),
                )
            }));
            context.queue = pq;
        }

        let (mod_vs, mod_fs) = split_long_edges_with_queue(
            &mut context.mesh,
            high,
            &adaptive_target_length,
            &mut feat_verts,
            &mut feat_edges,
            &mut context.queue,
            parameters.target_additional_vertices,
            true,
        );

        let new_vert_count = mod_vs.len();

        context.modified_vertices.extend(mod_vs);
        context.modified_faces.extend(mod_fs);

        equalize_valences(
            &mut context.mesh,
            parameters.minimum_triangle_angle,
            &feat_edges,
        );

        tangential_relaxation(&mut context.mesh, &_corners);

        if new_vert_count == 0 {
            break;
        }

        // TODO: implement reprojection.
    }

    context.mesh.vert_count() - input_verts
}

/// Split edges that exceed a given dynamic length. as computed by
/// `adaptive_target_length`.
#[allow(clippy::too_many_arguments)]
fn split_long_edges_with_queue<S, L, R>(
    mesh: &mut Redge<R>,
    high: S,
    adaptive_target_length: &L,
    feature_verts: &mut HashSet<VertId>,
    feature_edges: &mut HashSet<EdgeId>,
    pq: &mut PQueue<EdgeId, S>,
    iter_cap: usize,
    return_modified: bool,
) -> (Vec<VertId>, Vec<FaceId>)
where
    R: RedgeContainers,
    VertData<R>: InnerSpace<S>,
    S: RealField + num::traits::float::FloatCore + Mul<VertData<R>, Output = VertData<R>>,
    L: Fn(usize, usize) -> S,
{
    let mut mod_verts = Vec::new();
    let mut mod_faces = Vec::new();

    let mut iter_count = 0;
    while !pq.is_empty() && (iter_count < iter_cap || iter_cap == 0) {
        let (e_id, weight) = pq.pop().unwrap();
        let edge_handle = mesh.edge_handle(e_id);

        let adaptive_length =
            adaptive_target_length(edge_handle.v1().id().0, edge_handle.v2().id().0);
        let target_length = high * adaptive_length;
        debug_assert!(target_length > S::from(f32::EPSILON * 10.0).unwrap());
        if S::from(weight).unwrap() > target_length {
            // TODO(low): properly handle boundary edges in remeshing.
            if edge_handle.is_boundary() {
                continue;
            }

            let in_feature = feature_edges.contains(&edge_handle.id());

            let new_vert = mesh.split_edge(edge_handle.id());
            debug_assert!(correctness_state(mesh) == RedgeCorrectness::Correct);

            if return_modified {
                mod_verts.push(new_vert);

                mod_faces.extend(mesh.vert_handle(new_vert).incident_faces().map(|f| f.id()));
            }

            if in_feature {
                feature_verts.insert(new_vert);

                feature_edges.extend(mesh.vert_handle(new_vert).star_edges().map(|e| e.id()));
            }

            iter_count += 1;
        }
    }

    (mod_verts, mod_faces)
}

fn detect_features<S, R>(
    mesh: &mut Redge<R>,
    dihedral_angle_tolerance: S,
) -> (
    HashSet<VertId>, // Feature verts.
    HashSet<EdgeId>, // Feature edges.
    HashSet<VertId>, // Corner verts.
)
where
    R: RedgeContainers,
    S: RealField,
    VertData<R>: InnerSpace<S>,
{
    let is_feature_edge = |e: &EdgeHandle<R>| {
        let [o1, o2] = e.opposite_vertex_positions();
        let [v1, v2] = e.vertex_positions();

        let n1 = (v1.clone() - o1.clone())
            .cross(&(v2.clone() - o1))
            .normalized();

        let n2 = (v1 - o2.clone()).cross(&(v2 - o2)).normalized();

        let dihedral_angle = n1.dot(&n2).acos();

        S::from(dihedral_angle).unwrap() < dihedral_angle_tolerance
    };

    let mut feature_verts = HashSet::new();
    let mut feature_edges = HashSet::new();
    for edge in mesh.meta_edges() {
        if is_feature_edge(&edge) {
            feature_verts.insert(edge.v1().id());
            feature_verts.insert(edge.v2().id());

            feature_edges.insert(edge.id());
        }
    }

    // Find corners in the mesh.
    let mut corner_verts = HashSet::new();
    {
        let is_feature_edge = |e: &EdgeHandle<R>| feature_edges.contains(&(e.id()));

        for v in mesh.meta_verts() {
            let mut feature_edge_count = 0;
            for e in v.star_edges() {
                feature_edge_count += is_feature_edge(&e) as usize;
            }

            if feature_edge_count != 2 && feature_edge_count != 0 {
                corner_verts.insert(v.id());
            }
        }
    }

    (feature_verts, feature_edges, corner_verts)
}

fn equalize_valences<R, S>(mesh: &mut Redge<R>, minimum_angle: S, feature_edges: &HashSet<EdgeId>)
where
    R: RedgeContainers,
    VertData<R>: InnerSpace<S>,
    S: RealField + Mul<VertData<R>, Output = VertData<R>>,
{
    let is_feature_edge = |e: &EdgeHandle<R>| feature_edges.contains(&e.id());

    let target_val = |v: &VertHandle<R>| {
        if v.is_in_boundary() {
            4
        } else {
            6
        }
    };

    let deviation = |edge: &EdgeHandle<R>| {
        let mut res = 0;
        let a_valence = edge.v1().star_edges().count();
        let b_valence = edge.v2().star_edges().count();

        res += (a_valence as isize - target_val(&edge.v1()) as isize).abs()
            + (b_valence as isize - target_val(&edge.v2()) as isize).abs();

        for hedge in edge.hedge().radial_loop() {
            let valence = hedge.face_prev().source().star_edges().count();
            res += (valence as isize - target_val(&hedge.face_prev().source()) as isize).abs();
        }

        res
    };

    for i in 0..mesh.edges_meta.len() {
        if !mesh.edges_meta[i].is_active {
            continue;
        }

        let eid = EdgeId(i);
        let e = mesh.edge_handle(eid);
        // Don't flip if this would either break topology or flip a feature edge.
        if !e.can_flip() || e.is_boundary() || is_feature_edge(&e) {
            continue;
        }

        fn smallest_angle<S: RealField, V: VectorSpace<Scalar = S> + InnerSpace<S>>(
            v1: V,
            v2: V,
            v3: V,
        ) -> S {
            let d1 = (v2.clone() - v1.clone()).normalized();
            let d2 = (v3.clone() - v2.clone()).normalized();
            let d3 = (v1.clone() - v3.clone()).normalized();

            let a1 = (-d1.dot(&d2)).acos();
            let a2 = (-d2.dot(&d3)).acos();
            let a3 = (-d3.dot(&d1)).acos();

            a1.min(a2).min(a3)
        }

        let [o1, o2] = e.opposite_vertex_positions();
        let [v1, v2] = [e.v1().data().clone(), e.v2().data().clone()];
        let post_min_angle = smallest_angle(o1.clone(), o2.clone(), v1.clone())
            .min(smallest_angle(o1, o2, v2.clone()));

        if S::from(post_min_angle).unwrap() < minimum_angle {
            continue;
        }

        let deviation_pre = deviation(&e);

        mesh.flip_edge(eid);
        debug_assert!(correctness_state(mesh) == RedgeCorrectness::Correct);

        let e = mesh.edge_handle(eid);
        let deviation_post = deviation(&e);

        if deviation_pre <= deviation_post {
            mesh.flip_edge(eid);
        }
        debug_assert!(correctness_state(mesh) == RedgeCorrectness::Correct);
    }
}

fn tangential_relaxation<R, S>(mesh: &mut Redge<R>, corner_verts: &HashSet<VertId>)
where
    R: RedgeContainers,
    VertData<R>: InnerSpace<S>,
    S: RealField + Mul<VertData<R>, Output = VertData<R>>,
{
    let normals: Vec<_> = mesh
        .meta_faces()
        .map(|f| {
            let vs: Vec<_> = f.vertices().take(3).map(|v| v.data().clone()).collect();
            let d1 = vs[1].clone() - vs[0].clone();
            let d2 = vs[2].clone() - vs[0].clone();

            d1.cross(&d2).normalized()
        })
        .collect();

    let get_normal = |id: usize| normals[id].clone();

    let mut barycenters = vec![VertData::<R>::default(); mesh.vert_count()];
    let mut valences = vec![0; mesh.vert_count()];
    for vert in mesh.meta_verts() {
        let mut valence = 0;
        for neighbour in vert.neighbours() {
            barycenters[neighbour.id().0] += neighbour.data().clone();
            valence += 1;
        }
        valences[vert.id().0] = valence;
    }

    for (i, &valence) in valences.iter().enumerate() {
        barycenters[i] =
            barycenters[i].clone() * (S::from(1.0).unwrap() / S::from(valence).unwrap());
    }

    // This is way more readable the way it is.
    #[allow(clippy::needless_range_loop)]
    for i in 0..mesh.vert_data.len() {
        let vert = mesh.vert_data.get_mut(i as u64);
        if corner_verts.contains(&VertId(i)) {
            continue;
        }
        let normal = get_normal(i);
        let mut pos = vert.clone();
        pos = pos.clone() + normal.clone() * normal.dot(&(pos.clone() - barycenters[i].clone()));

        *vert = pos;
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Vector3;

    use crate::wavefront_loader::ObjData;

    use super::{incremental_refinement_with_context, *};

    // This is disabled as a test because it will chug our CI tools. But it is left here in case
    // verification is still needed.
    // #[test]
    fn _test_incremental_refinement_with_context() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/loop_cube.obj");

        // Convert to a type that admits arithmetic transformations.
        let vertices: Vec<_> = vertices
            .iter()
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

        let (vs, fs, _) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/before_refinement.obj");

        let number_of_subdivisions = 4000;
        let mut context = RemeshingContext::from_mesh(redge);
        for _i in 0..4 {
            incremental_refinement_with_context(
                &mut context,
                RemeshingParametersWithoutCollapse {
                    target_additional_vertices: number_of_subdivisions,
                    ..Default::default()
                },
                |_i, _j| 0.01,
                |&v| v,
            );
        }

        let (vs, fs, _) = context.mesh.to_face_list();
        ObjData::export(&(&vs, &fs), "out/after_refinement.obj");
    }
}
