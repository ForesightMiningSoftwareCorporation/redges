//! Topology handle around a half edge.
use crate::container_trait::RedgeContainers;
use crate::edge_handle::EdgeHandle;
use crate::face_handle::FaceHandle;
use crate::iterators::{FaceLoopHedgeIter, RadialHedgeIter};
use crate::vert_handle::VertHandle;

use crate::{HedgeId, HedgeMetaData, Redge};

/// Topology handle for a face.
pub struct HedgeHandle<'r, R: RedgeContainers> {
    id: HedgeId,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> HedgeHandle<'r, R> {
    pub(crate) fn new(id: HedgeId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    /// Id.
    pub fn id(&self) -> HedgeId {
        self.id
    }
    /// Vertex at the start of the half edge.
    pub fn source(&self) -> VertHandle<'r, R> {
        VertHandle::new(
            self.redge.hedges_meta[self.id.to_index()].source_id,
            self.redge,
        )
    }

    /// Vertex at the start of the end of the half edge.
    pub fn dest(&self) -> VertHandle<'r, R> {
        self.face_next().source()
    }

    /// Edge this half edge is incident on.
    pub fn edge(&self) -> EdgeHandle<'r, R> {
        EdgeHandle::new(self.metadata().edge_id, self.redge)
    }

    /// Next half edge in the radial cycle orbiting the edge.
    pub fn radial_next(&self) -> Self {
        debug_assert!(
            !self.metadata().radial_next_id.is_absent(),
            "Hedge {:?} radial next's is absent.",
            self.id
        );
        HedgeHandle::new(self.metadata().radial_next_id, self.redge)
    }
    /// Previous half edge in the radial cycle orbiting the edge.
    pub fn radial_prev(&self) -> Self {
        HedgeHandle::new(self.metadata().radial_prev_id, self.redge)
    }
    /// Next half edge in the face cycle.
    pub fn face_next(&self) -> Self {
        HedgeHandle::new(self.metadata().face_next_id, self.redge)
    }
    /// Previous half edge in the face cycle.
    pub fn face_prev(&self) -> Self {
        HedgeHandle::new(self.metadata().face_prev_id, self.redge)
    }
    /// Face this half edge belongs to.
    pub fn face(&self) -> FaceHandle<'r, R> {
        FaceHandle::new(self.metadata().face_id, self.redge)
    }
    /// Is this handle still pointing to a valid/active element?
    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &HedgeMetaData {
        &self.redge.hedges_meta[self.id.to_index()]
    }
    /// Iterator of the radial cycle of half edges orbiting the edge.
    pub fn radial_loop(&'r self) -> RadialHedgeIter<'r, R> {
        RadialHedgeIter::new(self.id(), self.redge)
    }
    /// Iterator of the cycle of half edges inside the face.
    pub fn face_loop(&'r self) -> FaceLoopHedgeIter<'r, R> {
        FaceLoopHedgeIter::new(self.id(), self.redge)
    }
}
