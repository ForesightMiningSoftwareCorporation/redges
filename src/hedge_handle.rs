use crate::*;
use edge_handle::*;
use face_handle::FaceHandle;
use vert_handle::*;

pub struct HedgeHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    id: HedgeId,
    redge: &'r mut Redge<V, E, F>,
}

impl<'r, V, E, F> HedgeHandle<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    pub(crate) fn new(id: HedgeId, redge: &'r mut Redge<V, E, F>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn id(&self) -> HedgeId {
        self.id
    }

    pub fn source(self) -> VertHandle<'r, V, E, F> {
        VertHandle::new(
            self.redge.hedges_meta[self.id.to_index()].source_id,
            self.redge,
        )
    }

    pub fn edge(self) -> EdgeHandle<'r, V, E, F> {
        EdgeHandle::new(self.my_ref().edge_id, self.redge)
    }

    pub fn radial_next(self) -> Self {
        HedgeHandle::new(self.my_ref().radial_next_id, self.redge)
    }

    pub fn radial_prev(self) -> Self {
        HedgeHandle::new(self.my_ref().radial_prev_id, self.redge)
    }

    pub fn face_next(self) -> Self {
        HedgeHandle::new(self.my_ref().face_next_id, self.redge)
    }

    pub fn face_prev(self) -> Self {
        HedgeHandle::new(self.my_ref().face_prev_id, self.redge)
    }

    pub fn face(self) -> FaceHandle<'r, V, E, F> {
        FaceHandle::new(self.my_ref().face_id, self.redge)
    }

    fn my_ref(&self) -> &HedgeMetaData {
        &self.redge.hedges_meta[self.id.to_index()]
    }
}
