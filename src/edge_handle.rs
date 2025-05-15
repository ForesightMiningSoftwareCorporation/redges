//! Topology handle around an edge.

use std::collections::HashSet;

use crate::container_trait::{EdgeData, RedgeContainers, VertData};
use crate::hedge_handle::HedgeHandle;
use crate::vert_handle::VertHandle;

use crate::HedgeId;
use crate::{
    container_trait::PrimitiveContainer, EdgeId, EdgeMetaData, Redge, StarCycleNode, VertId,
};

/// Topology handle for an edge.
pub struct EdgeHandle<'r, R: RedgeContainers> {
    id: EdgeId,
    pub(crate) redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> EdgeHandle<'r, R> {
    pub(crate) fn new(id: EdgeId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    /// Ids of the two endpoints.
    pub fn vertex_ids(&self) -> [VertId; 2] {
        self.metadata().vert_ids
    }

    /// positions of the two vertices.
    pub fn vertex_positions(&self) -> [VertData<R>; 2] {
        [self.v1().data().clone(), self.v2().data().clone()]
    }

    /// Id.
    pub fn id(&self) -> EdgeId {
        self.id
    }

    /// Handle to the first vertex of the edge.
    pub fn v1(&self) -> VertHandle<'r, R> {
        VertHandle::new(self.metadata().vert_ids[0], self.redge)
    }

    /// Handle to the second vertex of the edge.
    pub fn v2(&self) -> VertHandle<'r, R> {
        VertHandle::new(self.metadata().vert_ids[1], self.redge)
    }

    /// Handle to *a* half edge incident ont hsi edge.
    pub fn hedge(&self) -> HedgeHandle<'r, R> {
        HedgeHandle::new(self.metadata().hedge_id, self.redge)
    }

    /// Underlying data.
    pub fn data(&self) -> &EdgeData<R> {
        self.redge.edge_data.get(self.id.to_index() as u64)
    }

    /// Does this handle still point to a valid/active edge?
    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    pub(crate) fn metadata(&self) -> &EdgeMetaData {
        &self.redge.edges_meta[self.id.to_index()]
    }

    /// Identify which endpoint, if any, the given id corresponds to.
    pub fn vertex_endpoint(&self, vert_id: VertId) -> EdgeVertexType {
        let [v1, v2] = self.metadata().vert_ids;
        let v1_idx = v1.to_index();
        let v2_idx = v2.to_index();

        let index = vert_id.to_index();
        if index == v1_idx {
            EdgeVertexType::V1
        } else if index == v2_idx {
            EdgeVertexType::V2
        } else {
            EdgeVertexType::NotInEdge
        }
    }
    /// If the id is in the edge, returns the other endpoints ID. Can panic if
    /// the input is not either of the endpoints.
    pub fn opposite(&self, vert_id: VertId) -> VertId {
        match self.vertex_endpoint(vert_id) {
            EdgeVertexType::V1 => self.metadata().vert_ids[1],
            EdgeVertexType::V2 => self.metadata().vert_ids[0],
            EdgeVertexType::NotInEdge => panic!("Vertex id not part of this edge"),
        }
    }

    /// Get a handle to the iterator
    pub fn edge_cycle_at(&self, vert_id: VertId) -> StarCycleNode {
        match self.vertex_endpoint(vert_id) {
            EdgeVertexType::V1 => self.metadata().v1_cycle.clone(),
            EdgeVertexType::V2 => self.metadata().v2_cycle.clone(),
            EdgeVertexType::NotInEdge => panic!(
                "Vertex id {:?} not part of edge with id {:?}. That edge points to {:?}.",
                vert_id,
                self.id,
                self.vertex_ids(),
            ),
        }
    }

    /// Tests if collapsing the current edge would violate topology. If it
    /// returns true, collapsing the edge is safe.
    // TODO: this is assuming triangular faces right now.
    pub fn can_collapse(&self) -> bool {
        let v1 = self.v1();
        let v2 = self.v2();

        // If there's more than two vertices in the intersection of the neighbours of the endpoints of the edge,
        // then collapsing this edge will break topology. This happens when there is an edge connecting
        // the opposite vertices.
        // TODO: this logic was designed for triangular faces, it's not certain it works on polygonal ones.
        let set1: Vec<VertId> = v1.neighbours().map(|v| v.id()).collect();
        let count = v2.neighbours().filter(|v| set1.contains(&v.id())).count();

        if count > 2 {
            return false;
        }

        // Imagine, for example a triangulated strip:
        // *---*
        // |\  |
        // | \ |
        // |  \|
        // *---* <---- Bad idea to collapse this edge.
        // |\  |
        // | \ |
        // |  \|
        // *---*
        // Collapsing the middle edge would create degenerate geometry.
        let both_endpoints_are_in_boundary =
            self.v1().is_in_boundary() && self.v2().is_in_boundary();

        self.is_boundary() || !both_endpoints_are_in_boundary
    }

    /// Would flipping this edge break topology?
    pub fn can_flip(&self) -> bool {
        let [_, [t1, t2]] = self.butterfly_vertices_ids();
        let v1 = self.redge.vert_handle(t1);
        let v2 = self.redge.vert_handle(t2);

        let v1_set: HashSet<VertId> = HashSet::from_iter(v1.neighbours().map(|v| v.id()));
        let v2_set: HashSet<VertId> = HashSet::from_iter(v2.neighbours().map(|v| v.id()));

        !v1_set.contains(&v2.id()) && !v2_set.contains(&v1.id())
    }

    /// Is this edge in the boundary (has one or less incident faces).
    pub fn is_boundary(&self) -> bool {
        self.redge.edges_meta[self.id.to_index()].hedge_id == HedgeId::ABSENT
            || self.hedge().radial_loop().count() <= 1
    }

    /// Is there a half edge pointing to this edge?
    pub fn has_hedge(&self) -> bool {
        self.metadata().hedge_id != HedgeId::ABSENT
    }

    /// Only works on manifold geometry, and only on non-boundary edges.
    /// It returns the two vertices on the edge, followed by the two vertices
    /// transversal to it.
    pub(crate) fn butterfly_vertices_ids(&self) -> [[VertId; 2]; 2] {
        let [v1, v2] = self.vertex_ids();
        let v3 = self.hedge().face_prev().source().id();
        let v4 = self.hedge().radial_next().face_prev().source().id();

        [[v1, v2], [v3, v4]]
    }

    /// Get the positions of the two transverasl vertices for this edge (for triangular incident faces).
    pub(crate) fn opposite_vertex_positions(&self) -> [VertData<R>; 2] {
        let o1 = self.hedge().face_prev().source().data().clone();
        let o2 = self
            .hedge()
            .radial_next()
            .face_prev()
            .source()
            .data()
            .clone();

        [o1, o2]
    }
}

/// Relationship between the edge and a vertex id.
pub enum EdgeVertexType {
    /// First point in the edge.
    V1,
    /// Second point in the edge.
    V2,
    /// Point is not in the edge.
    NotInEdge,
}
