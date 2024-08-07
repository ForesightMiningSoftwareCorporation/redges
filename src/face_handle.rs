use container_trait::{FaceData, RedgeContainers};

use crate::hedge_handle::HedgeHandle;
use crate::vert_handle::VertHandle;

use crate::{container_trait::PrimitiveContainer, EdgeId, EdgeMetaData, Redge, VertId};

use crate::*;

pub struct FaceHandle<'r, R: RedgeContainers> {
    id: FaceId,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> FaceHandle<'r, R> {
    pub(crate) fn new(id: FaceId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn id(&self) -> FaceId {
        self.id
    }

    pub fn hedge(&self) -> HedgeHandle<'r, R> {
        HedgeHandle::new(self.metadata().hedge_id, self.redge)
    }

    pub fn data(&self) -> &FaceData<R::FaceContainer> {
        self.redge.face_data.get(self.id.to_index() as u64)
    }

    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &FaceMetaData {
        &self.redge.faces_meta[self.id.to_index()]
    }
}
