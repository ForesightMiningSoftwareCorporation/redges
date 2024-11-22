use crate::container_trait::RedgeContainers;
use crate::edge_handle::EdgeHandle;
use crate::face_handle::FaceHandle;
use crate::iterators::{FaceLoopHedgeIter, RadialHedgeIter};
use crate::vert_handle::VertHandle;

use crate::{HedgeId, HedgeMetaData, Redge};

pub struct HedgeHandle<'r, R: RedgeContainers> {
    id: HedgeId,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> HedgeHandle<'r, R> {
    pub(crate) fn new(id: HedgeId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn id(&self) -> HedgeId {
        self.id
    }

    pub fn source(&self) -> VertHandle<'r, R> {
        VertHandle::new(
            self.redge.hedges_meta[self.id.to_index()].source_id,
            self.redge,
        )
    }

    pub fn dest(&self) -> VertHandle<'r, R> {
        self.face_next().source()
    }

    pub fn edge(&self) -> EdgeHandle<'r, R> {
        EdgeHandle::new(self.metadata().edge_id, self.redge)
    }

    pub fn radial_next(&self) -> Self {
        debug_assert!(
            !self.metadata().radial_next_id.is_absent(),
            "Hedge {:?} radial next's is absent.",
            self.id
        );
        HedgeHandle::new(self.metadata().radial_next_id, self.redge)
    }

    pub fn radial_prev(&self) -> Self {
        HedgeHandle::new(self.metadata().radial_prev_id, self.redge)
    }

    pub fn face_next(&self) -> Self {
        HedgeHandle::new(self.metadata().face_next_id, self.redge)
    }

    pub fn face_prev(&self) -> Self {
        HedgeHandle::new(self.metadata().face_prev_id, self.redge)
    }

    pub fn face(&self) -> FaceHandle<'r, R> {
        FaceHandle::new(self.metadata().face_id, self.redge)
    }

    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &HedgeMetaData {
        &self.redge.hedges_meta[self.id.to_index()]
    }

    pub fn radial_loop(&'r self) -> RadialHedgeIter<'r, R> {
        RadialHedgeIter::new(self.id(), self.redge)
    }

    pub fn face_loop(&'r self) -> FaceLoopHedgeIter<'r, R> {
        FaceLoopHedgeIter::new(self.id(), self.redge)
    }
}
