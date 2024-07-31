use edge_handle::EdgeHandle;
use hedge_handle::HedgeHandle;
use vert_handle::VertHandle;

use crate::*;

pub struct FaceHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    id: FaceId,
    redge: &'r mut Redge<V, E, F>,
}

impl<'r, V, E, F> FaceHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    pub(crate) fn new(id: FaceId, redge: &'r mut Redge<V, E, F>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn id(&self) -> FaceId {
        self.id
    }

    pub fn hedge(self) -> HedgeHandle<'r, V, E, F> {
        HedgeHandle::new(self.my_ref().hedge_id, self.redge)
    }

    pub fn data(&self) -> &F::PrimitiveData {
        self.redge.face_data.get(self.id.to_index() as u64)
    }

    pub fn data_mut(&mut self) -> &mut F::PrimitiveData {
        self.redge.face_data.get_mut(self.id.to_index() as u64)
    }

    fn my_ref(&self) -> &FaceMetaData {
        &self.redge.faces_meta[self.id.to_index()]
    }
}
