//! This module implements quadric error based mesh simplification. It includes both
//! a variant that is agnostic to attributes and a version that measures attribute errors.

use core::f32;
use core::hash::Hash;
use std::{collections::BTreeSet, ops::Mul};

use crate::{
    algorithms::pqueue::PQueue,
    container_trait::{
        FaceAttributeGetter, FaceData, PrimitiveContainer, VertData, VertexAttributeGetter,
    },
    edge_handle::EdgeHandle,
    face_handle::{FaceHandle, FaceMetrics},
    mesh_deleter::MeshDeleter,
    validation::{correctness_state, RedgeCorrectness},
    wedge::WedgeDS,
    EdgeId, FaceId, VertId,
};
use crate::{container_trait::RedgeContainers, Redge};
use linear_isomorphic::prelude::*;
use nalgebra::{ComplexField, DMatrix, DVector, Matrix3, Vector3};
use num_traits::{
    float::{FloatCore, TotalOrder},
    Bounded, Float, Signed,
};

const EDGE_WEIGHT_PENALTY: f32 = 10.0;

/// Which kind of simplification strategy to use.
#[derive(Debug, PartialEq, Eq)]
pub enum SimplificationStrategy {
    /// Stop simplifying when no valid edges exist for preserving topology (closed meshes remain closed).
    Conservative,
    /// Do whatever it takes to simplify geoemtry for as long as there is geometry to simplify (may introduce boundaries on
    /// otherwise closed meshes).
    Aggressive,
}
/// Whether to simplify attributes as well as vertices.
#[derive(Debug, PartialEq, Eq)]
pub enum AttributeSimplification {
    /// Ignore attributes.
    NoAttributeSimplification,
    /// Use attributes in simplification.
    SimplifyAttributes,
}
/// What strategy to use to decide to stop simplificaiton.
pub enum StopCondition {
    /// Stop when a specified face cound is reached.
    FaceCount(usize),
    /// Stop when the error fals below this value.
    ErrorThreshold(f32),
}

/// Meta data to configure the simplification strategy.
pub struct QuadricSimplificationConfig {
    /// Whether to simplify aggressively or conservatively.
    pub strategy: SimplificationStrategy,
    /// Whether to consider attributes during simplification.
    pub attribute_simplification: AttributeSimplification,
    /// Specify when the simplificaitn should stop.
    pub stop_condition: StopCondition,
}

impl Default for QuadricSimplificationConfig {
    fn default() -> Self {
        Self {
            strategy: SimplificationStrategy::Conservative,
            attribute_simplification: AttributeSimplification::NoAttributeSimplification,
            stop_condition: StopCondition::FaceCount(10_000),
        }
    }
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

    let (mut mesh, cost) = match config.attribute_simplification {
        AttributeSimplification::SimplifyAttributes => {
            simplify_with_attributes(mesh, config, locked_vertex)
        }
        AttributeSimplification::NoAttributeSimplification => {
            simplify_without_attributes(mesh, config, locked_vertex)
        }
    };

    for i in 0..mesh.vert_data.len() {
        let pos = mesh.vert_data.get(i as u64).clone();
        *mesh.vert_data.get_mut(i as u64) = pos * (S::from(1.0).unwrap() / scale);
    }

    (mesh, cost)
}

fn simplify_with_attributes<S, R, L>(
    mesh: Redge<R>,
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
    let mut wedges = construct_wedges(&mesh);
    let attribute_count = mesh.face_data.get(0).attribute_count();

    // Start the queue with each edge's cost.
    let mut queue = initialize_queue_with_attributes::<S, R>(&mesh, &wedges, attribute_count);

    let mut deleter = crate::mesh_deleter::MeshDeleter::start_deletion(mesh);

    let mut worst_cost = <S as Float>::min_value();
    while !queue.is_empty() {
        let (
            QueueEdgeAttributeData {
                id: eid,
                geometric_cost,
                attributes: optimum,
                wedges_to_merge,
                wedge_order,
            },
            cost,
        ) = queue.pop().unwrap();

        // Check if the simplification has reached its goal.
        if match config.stop_condition {
            StopCondition::FaceCount(target_face_count) => {
                deleter.active_face_count() <= target_face_count
            }
            StopCondition::ErrorThreshold(threshold) => cost <= S::from(threshold).unwrap(),
        } {
            break;
        }

        let edge_handle = deleter.mesh().edge_handle(eid);

        // Skip inactive edges.
        if !edge_handle.is_active() {
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
        worst_cost = Float::max(worst_cost, geometric_cost);

        let edge_handle = deleter.mesh().edge_handle(eid);
        // Remove all edges that touch the current edge from the queue.
        for e in edge_handle
            .v1()
            .star_edges()
            .chain(edge_handle.v1().link_edges())
            .chain(edge_handle.v2().star_edges())
            .chain(edge_handle.v2().link_edges())
        {
            queue.remove(QueueEdgeAttributeData {
                id: e.id(),
                geometric_cost: S::default(),
                attributes: DVector::default(),
                wedge_order: Vec::new(),
                wedges_to_merge: Vec::new(),
            });
        }

        // Compute the edge collapse and update the new position (and attributes if needed).
        let edge_handle = deleter.mesh.edge_handle(eid);

        wedges.collapse_wedge(
            &optimum,
            &edge_handle,
            &wedges_to_merge,
            &wedge_order,
            attribute_count,
        );
        let v1 = edge_handle.v1().id();
        let v2 = edge_handle.v2().id();

        let vid = deleter.collapse_edge_and_fix(eid);

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
            if let Some(mismatched) = vids.iter().position(|id| !updated.contains(id)) {
                let data = deleter.mesh.face_data.get_mut(fid.to_index() as u64);
                data.attribute_vertices_mut()[mismatched] = vid;
            }
        }

        let vn = deleter.mesh.vert_handle(vid);
        let edges: Vec<_> = vn
            .star_edges()
            .chain(vn.link_edges())
            .map(|e| e.id())
            .collect();

        // The vertex could have no good neighbourhood. Consider for example
        // collapsing an edge in an isolated triangle. In these cases,
        // remove the geometry and move on.
        match edges.len() {
            0 => {
                deleter.remove_vert(vid);
                continue;
            }
            1 => {
                deleter.remove_edge(edges[0]);
                continue;
            }
            _ => {}
        }

        for eid in edges {
            let e = deleter.mesh.edge_handle(eid);
            assert!(e.is_active());

            let (mut cost, geom_cost, optimum, wedges_to_merge, wedge_order) =
                edge_cost_with_wedges(&e, &wedges, attribute_count);

            let [v1, v2] = e.vertex_ids();
            if locked_vertex(v1, &deleter.mesh) || locked_vertex(v2, &deleter.mesh) {
                cost += S::from(EDGE_WEIGHT_PENALTY * 100.).unwrap();
            }

            queue.push(
                QueueEdgeAttributeData {
                    id: e.id(),
                    geometric_cost: geom_cost,
                    attributes: optimum,
                    wedges_to_merge,
                    wedge_order,
                },
                cost,
            );
        }
    }

    if let StopCondition::FaceCount(target_face_count) = config.stop_condition {
        // Doing edge collapse after a certain point is very challenging, as a compromise,
        // if we reach here, and if we need a smaller mesh, we will just delete faces, if anyone wants
        // to try making an edge collapse that works no matter the situation you have my blessing.
        while deleter.active_face_count() > target_face_count
            && config.strategy == SimplificationStrategy::Aggressive
        {
            let face = deleter.mesh.faces_meta.iter().find(|f| f.is_active);
            if let Some(f) = face {
                deleter.remove_face(f.id);
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

fn simplify_without_attributes<S, R, L>(
    mesh: Redge<R>,
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
    // Start the queue with each edge's cost.
    let mut queue = initialize_queue_without_attributes::<S, R>(&mesh);

    let mut deleter = crate::mesh_deleter::MeshDeleter::start_deletion(mesh);

    let mut worst_cost = <S as Float>::min_value();
    while !queue.is_empty() {
        let (
            QueueEdgeSimpleData {
                id: eid,
                geometric_cost,
                optimum,
            },
            cost,
        ) = queue.pop().unwrap();

        // Check if the simplification has reached its goal.
        if match config.stop_condition {
            StopCondition::FaceCount(target_face_count) => {
                deleter.active_face_count() <= target_face_count
            }
            StopCondition::ErrorThreshold(threshold) => cost <= S::from(threshold).unwrap(),
        } {
            break;
        }

        println!("{}", cost);

        let edge_handle = deleter.mesh().edge_handle(eid);

        // Skip inactive edges.
        if !edge_handle.is_active() {
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

        debug_assert!(cost >= S::from(0.).unwrap(), "{}", cost);
        worst_cost = Float::max(worst_cost, geometric_cost);

        let edge_handle = deleter.mesh().edge_handle(eid);
        // Remove all edges that touch the current edge from the queue.
        for e in edge_handle
            .v1()
            .star_edges()
            .chain(edge_handle.v1().link_edges())
            .chain(edge_handle.v2().star_edges())
            .chain(edge_handle.v2().link_edges())
        {
            queue.remove(QueueEdgeSimpleData {
                id: e.id(),
                geometric_cost: S::default(),
                optimum: Vector3::default(),
            });
        }

        // Compute the edge collapse and update the new position (and attributes if needed).
        let vid = deleter.collapse_edge_and_fix(eid);

        deleter.mesh().vert_data(vid)[0] = optimum[0];
        deleter.mesh().vert_data(vid)[1] = optimum[1];
        deleter.mesh().vert_data(vid)[2] = optimum[2];

        let vn = deleter.mesh.vert_handle(vid);
        let edges: Vec<_> = vn
            .star_edges()
            .chain(vn.link_edges())
            .map(|e| e.id())
            .collect();

        // The vertex could have no good neighbourhood. Consider for example
        // collapsing an edge in an isolated triangle. In these cases,
        // remove the geometry and move on.
        match edges.len() {
            0 => {
                deleter.remove_vert(vid);
                continue;
            }
            1 => {
                deleter.remove_edge(edges[0]);
                continue;
            }
            _ => {}
        }

        for eid in edges {
            let e = deleter.mesh.edge_handle(eid);
            assert!(e.is_active());

            let (mut cost, geom_cost, optimum) = edge_cost_without_attributes(&e);
            debug_assert!(cost >= S::from(0.).unwrap());

            let [v1, v2] = e.vertex_ids();
            if locked_vertex(v1, &deleter.mesh) || locked_vertex(v2, &deleter.mesh) {
                cost += S::from(EDGE_WEIGHT_PENALTY * 100.).unwrap();
            }

            queue.push(
                QueueEdgeSimpleData {
                    id: e.id(),
                    geometric_cost: geom_cost,
                    optimum,
                },
                cost,
            );
        }
    }

    if let StopCondition::FaceCount(target_face_count) = config.stop_condition {
        // Doing edge collapse after a certain point is very challenging, as a compromise,
        // if we reach here, and if we need a smaller mesh, we will just delete faces, if anyone wants
        // to try making an edge collapse that works no matter the situation you have my blessing.
        while deleter.active_face_count() > target_face_count
            && config.strategy == SimplificationStrategy::Aggressive
        {
            let face = deleter.mesh.faces_meta.iter().find(|f| f.is_active);
            if let Some(f) = face {
                deleter.remove_face(f.id);
            }
        }
    }

    let (vert_frag, ..) = deleter.compute_fragmentation_maps();
    let mut res = deleter.end_deletion();
    debug_assert!(correctness_state(&res) == RedgeCorrectness::Correct);
    if res.face_data.get(0).attribute_count() >= 1 {
        // Update face connectivity upon termination.
        for face in 0..res.face_count() {
            let f = res.face_data.get_mut(face as u64);
            for v in f.attribute_vertices_mut() {
                *v = VertId(*vert_frag.get(&*v).unwrap());
            }
        }
    }

    (res, worst_cost)
}

fn face_area_after<R: RedgeContainers, S>(
    face: &FaceHandle<'_, R>,
    vert_id: VertId,
    p: &Vector3<S>,
) -> S
where
    VertData<R>: InnerSpace<S>,
    S: nalgebra::ComplexField + RealField,
{
    let mut verts = Vec::new();
    for vert in face.vertices() {
        let id = vert.id();
        let pos = if id == vert_id {
            *p
        } else {
            let d = vert.data().clone();
            Vector3::new(d[0], d[1], d[2])
        };

        verts.push(pos);
    }

    let d1 = verts[1] - verts[0];
    let d2 = verts[2] - verts[0];

    S::from_real(d1.cross(&d2).norm())
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

struct QueueEdgeAttributeData<S> {
    id: EdgeId,
    geometric_cost: S,
    attributes: DVector<S>,
    /// The wedges that must be merged after an edge collapse.
    wedges_to_merge: Vec<(usize, usize)>,
    /// The order in which the wedge attributes appear in the `attributes` vector.
    /// Needed when reconstructing face attribute values from a wedge.
    wedge_order: Vec<usize>,
}

impl<S> PartialEq for QueueEdgeAttributeData<S> {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl<S> Eq for QueueEdgeAttributeData<S> {}

impl<S> Hash for QueueEdgeAttributeData<S> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

struct QueueEdgeSimpleData<S> {
    id: EdgeId,
    geometric_cost: S,
    optimum: Vector3<S>,
}

impl<S> PartialEq for QueueEdgeSimpleData<S> {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl<S> Eq for QueueEdgeSimpleData<S> {}

impl<S> Hash for QueueEdgeSimpleData<S> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

/// For each vertex, compute the quadrics of its incident faces.
fn initialize_queue_with_attributes<S, R>(
    mesh: &Redge<R>,
    wedges: &WedgeDS<S>,
    attribute_count: usize,
) -> PQueue<QueueEdgeAttributeData<S>, S>
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
        let (cost, geom_cost, optimum, wedges_to_merge, wedge_order) =
            edge_cost_with_wedges(&edge, wedges, attribute_count);
        queue.push(
            QueueEdgeAttributeData {
                id: edge.id(),
                geometric_cost: geom_cost,
                attributes: optimum,
                wedges_to_merge,
                wedge_order,
            },
            cost,
        );
    }

    queue
}

fn initialize_queue_without_attributes<S, R>(mesh: &Redge<R>) -> PQueue<QueueEdgeSimpleData<S>, S>
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
        let (cost, geom_cost, optimum) = edge_cost_without_attributes(&edge);
        queue.push(
            QueueEdgeSimpleData {
                id: edge.id(),
                geometric_cost: geom_cost,
                optimum,
            },
            cost,
        );
    }

    queue
}

type EdgeCostData<S> = (S, S, nalgebra::DVector<S>, Vec<(usize, usize)>, Vec<usize>);
fn edge_cost_with_wedges<'r, S, R: RedgeContainers>(
    edge: &EdgeHandle<'r, R>,
    wedges: &WedgeDS<S>,
    attribute_count: usize,
) -> EdgeCostData<S>
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>> + ComplexField,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let (mut q, mut b, mut final_d, geom_q, geom_b, geom_d, wedges_to_merge, wedge_order) =
        wedges.wedge_quadric(edge, attribute_count);

    let mut add_phantom_plane =
        |e: &EdgeHandle<'r, R>, q: &mut DMatrix<S>, b: &mut DVector<S>, w: S| {
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

            top_left += nn * w;
            top_three += n * d * w;
            final_d += d * d * w;
        };

    // Add constraint planes for boundary edges.
    if edge.is_boundary() && edge.has_hedge() {
        // For each boundary edge that touches our active edge, we want
        // the final point to approximate the original boundary conditions as much as possible.
        // So add a phantom plane for each such edge, to nudge the final
        // edge position towards the boundary.
        let boundary_set: BTreeSet<_> = edge
            .v1()
            .star_edges()
            .chain(edge.v2().star_edges())
            .filter(|e| e.is_boundary() && e.has_hedge())
            .map(|e| e.id())
            .collect();
        for eid in boundary_set {
            add_phantom_plane(
                &edge.redge.edge_handle(eid),
                &mut q,
                &mut b,
                S::from(0.000001).unwrap(),
            );
        }
    }

    // When the edge is touching the boundary but it is not in it, we must be careful to
    // not move the boundary.
    let v1_is_fixed = edge.v1().is_in_boundary() && !edge.v2().is_in_boundary();
    let v2_is_fixed = edge.v2().is_in_boundary() && !edge.v1().is_in_boundary();

    let epsilon = 0.0001;
    let mut solve = |q: &mut DMatrix<S>, b: &mut DVector<S>| {
        let det = q.determinant();

        // Paolo Cignoni recommended that in cases where the quadrics are degenerate, we add
        // ghost planes which are orthogonal to the edges surrounding the collapsing edge.
        if det < S::from(epsilon).unwrap() {
            for e in edge.v1().star_edges().chain(edge.v2().star_edges()) {
                if !e.has_hedge() {
                    continue;
                }
                // TODO: maybe this should be proportional to the BB of the mesh or the average edge length.
                // instead of a cosntant value.
                add_phantom_plane(&e, q, b, S::from(0.00000001).unwrap());
            }
        }
        // Update the determinant.
        let det = q.determinant();

        // Safety check to prevent numerical problems when solving the system.
        if Float::abs(det) <= S::from(epsilon).unwrap() || v1_is_fixed || v2_is_fixed {
            return None;
        }
        q.clone().cholesky().map(|system| system.solve(&-b.clone()))
    };

    let (mut optimal, mut cost) = match solve(&mut q, &mut b) {
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

            let mut cost = (res.transpose() * &q * &res
                + b.transpose() * &res * S::from(2.).unwrap())[0]
                + final_d;

            if !cost.is_finite() {
                cost = S::from(0.).unwrap();
            }

            (
                res,
                cost + S::from(
                    (total_wedges as f32 + is_boundary as u8 as f32 * 10.) * EDGE_WEIGHT_PENALTY,
                )
                .unwrap(),
            )
        }
    };

    // If the optimum is very far away from the midpoint, it's likely that we had numerical issues in our matrix.
    // so set to the midpoint to be safe.
    let p = Vector3::new(optimal[0], optimal[1], optimal[2]);

    let mid = (edge.v1().data().clone() + edge.v2().data().clone()) * S::from(0.5).unwrap();
    let mid = Vector3::new(mid[0], mid[1], mid[2]);
    let v1 = edge.v1().data().clone();
    let v1 = Vector3::new(v1[0], v1[1], v1[2]);
    let v2 = edge.v2().data().clone();
    let v2 = Vector3::new(v2[0], v2[1], v2[2]);
    if (p - mid).norm() > (v1 - v2).norm() * S::from(4.).unwrap().real() {
        optimal[0] = mid[0];
        optimal[1] = mid[1];
        optimal[2] = mid[2];

        cost = (optimal.transpose() * &q * &optimal
            + b.transpose() * &optimal * S::from(2.).unwrap())[0]
            + final_d;

        if !cost.is_finite() {
            cost = S::from(0.).unwrap();
        }
    }

    // See if any of the faces becomes degenerate (i.e. area of 0).
    let f1 = if edge.has_hedge() {
        edge.hedge().face().id()
    } else {
        FaceId::ABSENT
    };
    let f2 = if edge.has_hedge() {
        edge.hedge().radial_next().face().id()
    } else {
        FaceId::ABSENT
    };
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

    let boundary_penalty = S::from(EDGE_WEIGHT_PENALTY).unwrap();
    assert!(cost.is_finite());
    let boundary_test =
        edge.is_boundary() || edge.v1().is_in_boundary() || edge.v2().is_in_boundary();
    let geom_cost = (optimal.fixed_rows::<3>(0).transpose() * geom_q * optimal.fixed_rows::<3>(0)
        + geom_b.transpose() * optimal.fixed_rows::<3>(0) * S::from(2.).unwrap())[0]
        + geom_d;
    (
        cost + S::from(boundary_test as u8 as f32).unwrap() * boundary_penalty,
        geom_cost,
        optimal,
        wedges_to_merge,
        wedge_order,
    )
}

fn edge_cost_without_attributes<'r, S, R: RedgeContainers>(
    edge: &EdgeHandle<'r, R>,
) -> (S, S, Vector3<S>)
where
    S: RealField + Mul<VertData<R>, Output = VertData<R>> + ComplexField,
    VertData<R>: InnerSpace<S> + VertexAttributeGetter<S>,
    FaceData<R>: FaceAttributeGetter<S>,
{
    let edge_quadric = |edge: &EdgeHandle<'r, R>| {
        // In non-manifold meshes, sometimes we generate weird edges that don't have incident faces.
        // but they could still connect to full faces whose quadrics we must optimize.
        let [f1, f2] = if edge.has_hedge() {
            let f1 = edge.hedge().face().id();
            let f2 = edge.hedge().radial_next().face().id();

            [f1, f2]
        } else {
            [FaceId::ABSENT, FaceId::ABSENT]
        };

        let mut final_a = Matrix3::<S>::zeros();
        let mut final_b = Vector3::<S>::zeros();
        let mut final_c = S::zero();
        for (_vert, face) in edge
            .v1()
            .incident_faces()
            .map(|f| (edge.v1(), f))
            .chain(edge.v2().incident_faces().map(|f| (edge.v2(), f)))
            // Skip the two faces sharing the edge.
            .filter(|f| f.1.id() != f1 && f.1.id() != f2)
        {
            let (a, b, c) = crate::wedge::face_geometric_quadric(&face);
            let area = face.area();

            final_a += a * area;
            final_b += b * area;
            final_c += c * area;
        }

        (final_a, final_b, final_c)
    };

    let (mut q, mut b, mut final_d) = edge_quadric(edge);

    let mut add_phantom_plane =
        |e: &EdgeHandle<'r, R>, q: &mut Matrix3<S>, b: &mut Vector3<S>, w: S| {
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

            *q += nn * w;
            *b += n * d * w;
            final_d += d * d * w;
        };

    // Add constraint planes for boundary edges.
    if edge.is_boundary() && edge.has_hedge() {
        // For each boundary edge that touches our active edge, we want
        // the final point to approximate the original boundary conditions as much as possible.
        // So add a phantom plane for each such edge, to nudge the final
        // edge position towards the boundary.
        let boundary_set: BTreeSet<_> = edge
            .v1()
            .star_edges()
            .chain(edge.v2().star_edges())
            .filter(|e| e.is_boundary() && e.has_hedge())
            .map(|e| e.id())
            .collect();
        for eid in boundary_set {
            add_phantom_plane(
                &edge.redge.edge_handle(eid),
                &mut q,
                &mut b,
                S::from(0.000001).unwrap(),
            );
        }
    }

    // When the edge is touching the boundary but it is not in it, we must be careful to
    // not move the boundary.
    let v1_is_fixed = edge.v1().is_in_boundary() && !edge.v2().is_in_boundary();
    let v2_is_fixed = edge.v2().is_in_boundary() && !edge.v1().is_in_boundary();

    let epsilon = 0.0001;
    let mut solve = |q: &mut Matrix3<S>, b: &mut Vector3<S>| {
        let det = q.determinant();

        // Paolo Cignoni recommended that in cases where the quadrics are degenerate, we add
        // ghost planes which are orthogonal to the edges surrounding the collapsing edge.
        if det < S::from(epsilon).unwrap() {
            for e in edge.v1().star_edges().chain(edge.v2().star_edges()) {
                if !e.has_hedge() {
                    continue;
                }
                // TODO: maybe this should be proportional to the BB of the mesh or the average edge length.
                // instead of a cosntant value.
                add_phantom_plane(&e, q, b, S::from(0.00000001).unwrap());
            }
        }
        // Update the determinant.
        let det = q.determinant();

        // Safety check to prevent numerical problems when solving the system.
        if Float::abs(det) <= S::from(epsilon).unwrap() || v1_is_fixed || v2_is_fixed {
            return None;
        }
        (*q).cholesky().map(|system| system.solve(&-*b))
    };

    let (mut optimal, mut cost) = match solve(&mut q, &mut b) {
        Some(optimal) => {
            let cost = (optimal.transpose() * q * optimal
                + b.transpose() * optimal * S::from(2.).unwrap())[0]
                + final_d;

            // TODO: is this caused just by floating point issues, or do we have a bug?
            let cost = cost.max(S::from(0.0).unwrap());

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
            let mut res = Vector3::zeros();

            // Penalize collapsing edges with larger numbers of incident wedges at their endpoints.
            let is_boundary = edge.v1().is_in_boundary() || edge.v2().is_in_boundary();

            res[0] = target[0];
            res[1] = target[1];
            res[2] = target[2];

            let cost = (res.transpose() * q * res + b.transpose() * res * S::from(2.).unwrap())[0]
                + final_d;
            // TODO: Is this just numerical problems or are the above values not wll made?
            let mut cost = S::from(0.).unwrap().max(cost);

            if !cost.is_finite() {
                cost = S::from(0.).unwrap();
            }

            (
                res,
                cost + S::from((is_boundary as u8 as f32 * 10.) * EDGE_WEIGHT_PENALTY).unwrap(),
            )
        }
    };

    // If the optimum is very far away from the midpoint, it's likely that we had numerical issues in our matrix.
    // so set to the midpoint to be safe.
    let p = Vector3::new(optimal[0], optimal[1], optimal[2]);

    let mid = (edge.v1().data().clone() + edge.v2().data().clone()) * S::from(0.5).unwrap();
    let mid = Vector3::new(mid[0], mid[1], mid[2]);
    let v1 = edge.v1().data().clone();
    let v1 = Vector3::new(v1[0], v1[1], v1[2]);
    let v2 = edge.v2().data().clone();
    let v2 = Vector3::new(v2[0], v2[1], v2[2]);
    if (p - mid).norm() > (v1 - v2).norm() * S::from(4.).unwrap().real() {
        optimal[0] = mid[0];
        optimal[1] = mid[1];
        optimal[2] = mid[2];

        cost = (optimal.transpose() * q * optimal + b.transpose() * optimal * S::from(2.).unwrap())
            [0]
            + final_d;
        // TODO: is this caused just by floating point issues, or do we have a bug?
        cost = cost.max(S::from(0.0).unwrap());

        if !cost.is_finite() {
            cost = S::from(0.).unwrap();
        }
    }

    // See if any of the faces becomes degenerate (i.e. area of 0).
    let f1 = if edge.has_hedge() {
        edge.hedge().face().id()
    } else {
        FaceId::ABSENT
    };
    let f2 = if edge.has_hedge() {
        edge.hedge().radial_next().face().id()
    } else {
        FaceId::ABSENT
    };
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

    let boundary_penalty = S::from(EDGE_WEIGHT_PENALTY).unwrap();
    debug_assert!(cost.is_finite());
    debug_assert!(cost >= S::from(0.).unwrap());

    let boundary_test =
        edge.is_boundary() || edge.v1().is_in_boundary() || edge.v2().is_in_boundary();
    let geom_cost = (optimal.fixed_rows::<3>(0).transpose() * q * optimal.fixed_rows::<3>(0)
        + b.transpose() * optimal.fixed_rows::<3>(0) * S::from(2.).unwrap())[0]
        + final_d;
    (
        cost + S::from(boundary_test as u8 as f32).unwrap() * boundary_penalty,
        geom_cost,
        optimal,
    )
}

/// Test if an edge collapse would flip the direction of a face normal.
fn collapse_would_flip_normal<R, S>(edge: &EdgeHandle<'_, R>, new_pos: &VertData<R>) -> bool
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
        let new_normal = get_normal(new_pos, &pp, &pn).normalized();

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

    use crate::wavefront_loader::ObjData;

    use super::*;

    // Let's not make it a test for now, since we would need to add the assets.
    #[test]
    fn _test_quadric_simplification_closed() {
        let ObjData {
            vertices,
            vertex_face_indices,
            ..
        } = ObjData::from_disk_file("assets/armadillo.obj");
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
        let target = redge.face_count() / 2;
        let (redge, _cost) = quadric_simplify(
            redge,
            QuadricSimplificationConfig {
                strategy: SimplificationStrategy::Aggressive,
                attribute_simplification: AttributeSimplification::NoAttributeSimplification,
                stop_condition: StopCondition::FaceCount(target),
            },
            |_, _| false,
        );

        let (vs, fs, _) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_closed.obj");
    }

    // Let's not make it a test for now, since we would need to add the assets.
    // #[test]
    fn _test_quadric_simplification_boundary() {
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

        let target = redge.face_count() / 10;
        let start = Instant::now();
        let (redge, _cost) = quadric_simplify(
            redge,
            QuadricSimplificationConfig {
                strategy: SimplificationStrategy::Conservative,
                attribute_simplification: AttributeSimplification::NoAttributeSimplification,
                stop_condition: StopCondition::FaceCount(target),
            },
            |_, _| false,
        );
        let duration = start.elapsed();
        println!("Time elapsed in simplification is: {:?}", duration);

        let (vs, fs, _) = redge.to_face_list();
        ObjData::export(&(&vs, &fs), "out/simplified_flat_donut.obj");
    }
}
