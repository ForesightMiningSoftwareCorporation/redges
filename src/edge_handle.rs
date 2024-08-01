use crate::hedge_handle::HedgeHandle;
use crate::vert_handle::VertHandle;

use crate::{EdgeId, EdgeMetaData, PrimitiveContainer, Redge, VertId};

pub struct EdgeHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    id: EdgeId,
    redge: &'r mut Redge<V, E, F>,
}

impl<'r, V, E, F> EdgeHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    pub(crate) fn new(id: EdgeId, redge: &'r mut Redge<V, E, F>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn vertex_ids(&self) -> [VertId; 2] {
        self.my_ref().vert_ids
    }

    pub fn id(&self) -> EdgeId {
        self.id
    }

    pub fn v0(self) -> VertHandle<'r, V, E, F> {
        VertHandle::new(self.my_ref().vert_ids[0], self.redge)
    }

    pub fn v1(self) -> VertHandle<'r, V, E, F> {
        VertHandle::new(self.my_ref().vert_ids[1], self.redge)
    }

    pub fn hedge(self) -> HedgeHandle<'r, V, E, F> {
        HedgeHandle::new(self.my_ref().hedge_id, self.redge)
    }

    pub fn data(&self) -> &E::PrimitiveData {
        self.redge.edge_data.get(self.id.to_index() as u64)
    }

    pub fn data_mut(&mut self) -> &mut E::PrimitiveData {
        self.redge.edge_data.get_mut(self.id.to_index() as u64)
    }

    fn my_ref(&self) -> &EdgeMetaData {
        &self.redge.edges_meta[self.id.to_index()]
    }
}
