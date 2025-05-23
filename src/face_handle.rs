//! Topology handle around a face.
use container_trait::{FaceData, RedgeContainers};
use iterators::{FaceLoopHedgeIter, FaceVertIterator};
use linear_isomorphic::prelude::*;
use std::collections::BTreeSet;
use std::fmt::Debug;

use crate::hedge_handle::HedgeHandle;
use crate::vert_handle::VertHandle;

use crate::{container_trait::PrimitiveContainer, Redge};

use crate::*;

/// Topology handle for a face.
pub struct FaceHandle<'r, R: RedgeContainers> {
    id: FaceId,
    pub(crate) redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> FaceHandle<'r, R> {
    pub(crate) fn new(id: FaceId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    /// Id.
    pub fn id(&self) -> FaceId {
        self.id
    }

    /// A random half edge within the face.
    pub fn hedge(&self) -> HedgeHandle<'r, R> {
        HedgeHandle::new(self.metadata().hedge_id, self.redge)
    }

    /// Underlying data.
    pub fn data(&self) -> &FaceData<R> {
        self.redge.face_data.get(self.id.to_index() as u64)
    }

    /// Is the face still valid/active.
    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &FaceMetaData {
        &self.redge.faces_meta[self.id.to_index()]
    }

    /// Vertices in the face.
    #[inline]
    pub fn vertices(&'r self) -> impl Iterator<Item = VertHandle<'r, R>> + 'r {
        FaceVertIterator {
            _redge: self.redge,
            face_loop: FaceLoopHedgeIter::new(self.hedge().id(), self.redge),
        }
    }

    /// Ids of the vertices in the face.
    #[inline]
    pub fn vertex_ids(&'r self) -> impl Iterator<Item = VertId> + 'r {
        self.vertices().map(|v| v.id())
    }

    /// If part of this face, return a handle to the vertex, otherwise `None`.
    #[inline]
    pub fn vertex_by_id(&'r self, vid: VertId) -> Option<VertHandle<'r, R>> {
        if self.hedge().face_loop().any(|h| h.source().id() == vid) {
            return Some(self.redge.vert_handle(vid));
        }

        None
    }

    /// If this face doesn't share an edge with any other one it will return true.
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

/// Which issues can be present in a face.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FaceDegeneracies {
    /// This face is topologically adequate.
    None,
    /// Two sides.
    Digon,
    /// One side.
    Monogon,
    /// Two triangles sharing more than two vertex indices.
    Doppelganger,
}

/// Trait allowing geometry queries in a face.
pub trait FaceMetrics<V, N, S>
where
    V: Clone + Debug,
    S: RealField,
{
    /// Face area.
    fn area(&self) -> S;
    /// Normal to the face, assumes a manifold mesh embedded in R^3.
    fn normal(&self) -> N;
    /// Unit normal to the face, assumes a manifold mesh embedded in R^3.
    fn unit_normal(&self) -> N;
    /// Average of the face vertex positions.
    fn centroid(&self) -> V;
}

impl<R: RedgeContainers, S> FaceMetrics<VertData<R>, VertData<R>, S> for FaceHandle<'_, R>
where
    VertData<R>: InnerSpace<S>,
    S: RealField,
{
    // TODO: only works on triangular faces for now.
    fn area(&self) -> S {
        self.normal().norm() / S::from(2.0).unwrap()
    }

    fn normal(&self) -> VertData<R> {
        // This implementation is expensive but reduces nuemrical errors by grabbing the best pair of edges to
        // compute a normal with.
        let pts: Vec<_> = self.vertices().map(|v| v.data().clone()).collect();
        let best_vertex = angle_closest_to_90(self);
        let d1 = pts[(best_vertex + 2) % 3].clone() - pts[(best_vertex + 1) % 3].clone();
        let d2 = pts[best_vertex].clone() - pts[(best_vertex + 1) % 3].clone();

        d1.cross(&d2)
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

/// Helper function to identify which vertex of a face has the best angle for numerics.
pub fn angle_closest_to_90<R: RedgeContainers, S>(face: &FaceHandle<'_, R>) -> usize
where
    VertData<R>: InnerSpace<S>,
    S: RealField,
{
    let points: Vec<_> = face.vertices().map(|v| v.data().clone()).collect();
    let mut best_cos = S::from(2.0).unwrap();
    let mut selected_i = 0;

    for i in 0..3 {
        let d1 = (points[i].clone() - points[(i + 1) % 3].clone()).normalized();
        let d2 = (points[(i + 2) % 3].clone() - points[(i + 1) % 3].clone()).normalized();

        let cos_abs = d1.dot(&d2).abs();

        if cos_abs < best_cos {
            best_cos = cos_abs;
            selected_i = i;
        }
    }

    selected_i
}
