use crate::container_trait::{EdgeData, RedgeContainers};
use crate::hedge_handle::HedgeHandle;
use crate::vert_handle::VertHandle;

use crate::{
    container_trait::PrimitiveContainer, EdgeId, EdgeMetaData, Redge, StarCycleNode, VertId,
};

pub struct EdgeHandle<'r, R: RedgeContainers> {
    id: EdgeId,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> EdgeHandle<'r, R> {
    pub(crate) fn new(id: EdgeId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn vertex_ids(&self) -> [VertId; 2] {
        self.metadata().vert_ids
    }

    pub fn id(&self) -> EdgeId {
        self.id
    }

    pub fn v1(&self) -> VertHandle<'r, R> {
        VertHandle::new(self.metadata().vert_ids[0], self.redge)
    }

    pub fn v2(&self) -> VertHandle<'r, R> {
        VertHandle::new(self.metadata().vert_ids[1], self.redge)
    }

    pub fn hedge(&self) -> HedgeHandle<'r, R> {
        HedgeHandle::new(self.metadata().hedge_id, self.redge)
    }

    pub fn data(&self) -> &EdgeData<R::EdgeContainer> {
        self.redge.edge_data.get(self.id.to_index() as u64)
    }

    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &EdgeMetaData {
        &self.redge.edges_meta[self.id.to_index()]
    }

    pub fn vertex_endpoint(&self, vert_id: VertId) -> EdgeVertexType {
        let [v1, v2] = self.metadata().vert_ids;
        let v1_idx = v1.to_index();
        let v2_idx = v2.to_index();

        match vert_id.to_index() {
            v1_idx => EdgeVertexType::V1,
            v2_idx => EdgeVertexType::V2,
            _ => EdgeVertexType::NotInEdge,
        }
    }

    pub fn opposite(&self, vert_id: VertId) -> VertId {
        match self.vertex_endpoint(vert_id) {
            EdgeVertexType::V1 => self.metadata().vert_ids[0],
            EdgeVertexType::V2 => self.metadata().vert_ids[1],
            EdgeVertexType::NotInEdge => panic!("Vertex id not part of this edge"),
        }
    }

    pub fn edge_cycle_at(&self, vert_id: VertId) -> StarCycleNode {
        match self.vertex_endpoint(vert_id) {
            EdgeVertexType::V1 => self.metadata().v1_cycle.clone(),
            EdgeVertexType::V2 => self.metadata().v2_cycle.clone(),
            EdgeVertexType::NotInEdge => panic!("Vertex id not part of this edge"),
        }
    }

    pub(crate) fn vertex_cycle_pointers_v1(&self) -> StarCycleNode {
        self.metadata().v1_cycle.clone()
    }

    pub(crate) fn vertex_cycle_pointers_v2(&self) -> StarCycleNode {
        self.metadata().v2_cycle.clone()
    }
}

pub enum EdgeVertexType {
    V1,
    V2,
    NotInEdge,
}
