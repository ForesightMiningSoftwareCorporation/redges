use container_trait::{FaceData, RedgeContainers};
use iterators::{FaceLoopHedgeIter, FaceVertIterator};
use linear_isomorphic::prelude::*;
use std::collections::{BTreeSet, HashSet};
use std::fmt::Debug;

use crate::hedge_handle::HedgeHandle;
use crate::vert_handle::VertHandle;

use crate::{container_trait::PrimitiveContainer, Redge};

use crate::*;

pub struct FaceHandle<'r, R: RedgeContainers> {
    id: FaceId,
    pub(crate) redge: &'r Redge<R>,
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

    pub fn data(&self) -> &FaceData<R> {
        self.redge.face_data.get(self.id.to_index() as u64)
    }

    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &FaceMetaData {
        &self.redge.faces_meta[self.id.to_index()]
    }

    #[inline]
    pub fn vertices(&'r self) -> impl Iterator<Item = VertHandle<'r, R>> + '_ {
        FaceVertIterator {
            _redge: self.redge,
            face_loop: FaceLoopHedgeIter::new(self.hedge().id(), self.redge),
        }
    }

    #[inline]
    pub fn vertex_ids(&'r self) -> impl Iterator<Item = VertId> + '_ {
        self.vertices().map(|v| v.id())
    }

    #[inline]
    pub fn vertex_by_id(&'r self, vid: VertId) -> Option<VertHandle<'r, R>> {
        if self.hedge().face_loop().any(|h| h.source().id() == vid) {
            return Some(self.redge.vert_handle(vid));
        }

        None
    }

    #[inline]
    pub fn is_isolated(&'r self) -> bool {
        self.hedge()
            .face_loop()
            .all(|h| h.radial_loop().count() == 1)
    }

    /// Count the number of sides in the face.
    pub fn side_count(&self) -> usize {
        self.hedge().face_loop().count()
    }

    /// Inspects the face for topological degeneracies such as a face with two sides or
    /// two faces sharing more than two vertices.
    pub fn check_degeneracies(&self) -> FaceDegeneracies {
        match self.side_count() {
            0 => panic!("Encountered a face with no sides."),
            1 => return FaceDegeneracies::Monogon,
            2 => return FaceDegeneracies::Digon,
            _ => {}
        }

        let s0: BTreeSet<_> = self.vertex_ids().collect();

        for h in self.hedge().face_loop() {
            for hedge in h.radial_loop().filter(|o| o.id() != h.id()) {
                let s1: BTreeSet<_> = hedge.face().vertex_ids().collect();

                let intersect = s0.intersection(&s1);
                if intersect.count() > 2 {
                    return FaceDegeneracies::Doppelganger;
                }
            }
        }

        FaceDegeneracies::None
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FaceDegeneracies {
    None,
    /// Two sides.
    Digon,
    /// One side.
    Monogon,
    /// Two triangles sharing more than two vertex indices.
    Doppelganger,
}

pub trait FaceMetrics<V, N, S>
where
    V: Clone + Debug,
    S: RealField,
{
    fn area(&self) -> S;
    fn normal(&self) -> N;
    fn unit_normal(&self) -> N;
    fn centroid(&self) -> V;
}

impl<'r, R: RedgeContainers, S> FaceMetrics<VertData<R>, VertData<R>, S> for FaceHandle<'r, R>
where
    VertData<R>: InnerSpace<S>,
    S: RealField,
{
    // TODO: only works on triangular faces for now.
    fn area(&self) -> S {
        self.normal().norm() / S::from(2.0).unwrap()
    }

    fn normal(&self) -> VertData<R> {
        let p1 = self.hedge().source().data().clone();
        let p2 = self.hedge().face_next().source().data().clone();
        let p3 = self.hedge().face_prev().source().data().clone();

        let e1 = p2 - p1.clone();
        let e2 = p3 - p1;

        e1.cross(&e2)
    }

    fn unit_normal(&self) -> VertData<R> {
        self.normal().normalized()
    }

    fn centroid(&self) -> VertData<R> {
        let mut centroid = VertData::<R>::default();
        let mut count = 0;
        for v in self.vertices() {
            centroid += v.data().clone();
            count += 1;
        }

        centroid * (S::from(1.0).unwrap() / S::from(count).unwrap())
    }
}
