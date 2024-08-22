use crate::container_trait::{RedgeContainers, VertData};
use crate::edge_handle::EdgeHandle;

use crate::iterators::{VertexStarEdgesIter, VertexStarVerticesIter};
use crate::{container_trait::PrimitiveContainer, Redge, VertId, VertMetaData};

pub struct VertHandle<'r, R: RedgeContainers> {
    id: VertId,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertHandle<'r, R> {
    pub(crate) fn new(id: VertId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn id(&self) -> VertId {
        self.id
    }

    pub fn edge(&self) -> EdgeHandle<'r, R> {
        EdgeHandle::new(
            self.redge.verts_meta[self.id.to_index()].edge_id,
            self.redge,
        )
    }

    pub fn data(&self) -> &VertData<R> {
        self.redge.vert_data.get(self.id.to_index() as u64)
    }

    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &VertMetaData {
        &self.redge.verts_meta[self.id.to_index()]
    }

    pub fn neighbours(&self) -> VertexStarVerticesIter<'r, R> {
        VertexStarVerticesIter::new(self.id, self.redge)
    }

    pub fn star_edges(&self) -> VertexStarEdgesIter<'r, R> {
        VertexStarEdgesIter::new(self.id, self.redge)
    }
}
