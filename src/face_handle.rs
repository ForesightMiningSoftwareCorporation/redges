use crate::hedge_handle::HedgeHandle;
use crate::vert_handle::VertHandle;

use crate::{EdgeId, EdgeMetaData, PrimitiveContainer, Redge, VertId};

use crate::*;

pub struct FaceHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    id: FaceId,
    redge: &'r Redge<V, E, F>,
}

impl<'r, V, E, F> FaceHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    pub(crate) fn new(id: FaceId, redge: &'r Redge<V, E, F>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn id(&self) -> FaceId {
        self.id
    }

    pub fn hedge(&self) -> HedgeHandle<'r, V, E, F> {
        HedgeHandle::new(self.metadata().hedge_id, self.redge)
    }

    pub fn data(&self) -> &F::PrimitiveData {
        self.redge.face_data.get(self.id.to_index() as u64)
    }

    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &FaceMetaData {
        &self.redge.faces_meta[self.id.to_index()]
    }
}
