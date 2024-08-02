use crate::edge_handle::EdgeHandle;
use crate::hedge_handle::HedgeHandle;

use crate::{EdgeId, PrimitiveContainer, Redge, VertId, VertMetaData};

pub struct VertHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    id: VertId,
    redge: &'r mut Redge<V, E, F>,
}

impl<'r, V, E, F> VertHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    pub(crate) fn new(id: VertId, redge: &'r mut Redge<V, E, F>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn id(&self) -> VertId {
        self.id
    }

    pub fn edge(self) -> EdgeHandle<'r, V, E, F> {
        EdgeHandle::new(
            self.redge.verts_meta[self.id.to_index()].edge_id,
            self.redge,
        )
    }

    pub fn data(&self) -> &V::PrimitiveData {
        self.redge.vert_data.get(self.id.to_index() as u64)
    }

    pub fn data_mut(&mut self) -> &mut V::PrimitiveData {
        self.redge.vert_data.get_mut(self.id.to_index() as u64)
    }

    pub fn is_active(&self) -> bool {
        self.my_ref().is_active
    }

    fn my_ref(&self) -> &VertMetaData {
        &self.redge.verts_meta[self.id.to_index()]
    }
}
