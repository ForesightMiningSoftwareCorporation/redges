use std::collections::{BTreeSet, HashSet};

use crate::{
    container_trait::RedgeContainers, edge_handle::EdgeHandle, face_handle::FaceHandle,
    hedge_handle::HedgeHandle, vert_handle::VertHandle, EdgeId, FaceId, HedgeId, Redge, VertId,
};

pub struct VertexStarVerticesIter<'r, R: RedgeContainers> {
    focused_vertex: VertId,
    iterator: VertexStarEdgesIter<'r, R>,

    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertexStarVerticesIter<'r, R> {
    pub(crate) fn new(vert_id: VertId, redge: &'r Redge<R>) -> Self {
        VertexStarVerticesIter {
            focused_vertex: vert_id,
            iterator: VertexStarEdgesIter::new(vert_id, redge),
            redge: redge,
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for VertexStarVerticesIter<'r, R> {
    type Item = VertHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        let res = self
            .iterator
            .next()
            .map(|e| self.redge.vert_handle(e.opposite(self.focused_vertex)));

        res
    }
}

pub struct VertexStarEdgesIter<'r, R: RedgeContainers> {
    pub(crate) start_edge: EdgeId,
    pub(crate) current_edge: EdgeId,
    pub(crate) focused_vertex: VertId,
    pub(crate) start: bool,

    pub(crate) redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertexStarEdgesIter<'r, R> {
    pub(crate) fn new(vert_id: VertId, redge: &'r Redge<R>) -> Self {
        VertexStarEdgesIter {
            focused_vertex: vert_id,
            start_edge: redge.verts_meta[vert_id.to_index()].edge_id,
            current_edge: redge.verts_meta[vert_id.to_index()].edge_id,
            redge: redge,
            start: true,
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for VertexStarEdgesIter<'r, R> {
    type Item = EdgeHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_edge == EdgeId::ABSENT {
            return None;
        }
        let edge = self.redge.edge_handle(self.current_edge);
        let edge_cycle = edge.edge_cycle_at(self.focused_vertex);

        if self.current_edge == self.start_edge {
            if !self.start {
                return None;
            } else {
                self.start = false;
            }
        }

        self.current_edge = edge_cycle.next_edge;

        Some(EdgeHandle::new(edge.id(), self.redge))
    }
}

// Note that the implementation of this iterator could be much faster if we assumed
// manifold topology.
pub struct VertexLinkEdgesIter<'r, R: RedgeContainers> {
    focused_vertex: VertId,
    edge_iter: VertexStarEdgesIter<'r, R>,
    orbit_iter: RadialHedgeIter<'r, R>,
    seen: BTreeSet<EdgeId>,
}

impl<'r, R: RedgeContainers> VertexLinkEdgesIter<'r, R> {
    pub(crate) fn new(vid: VertId, redge: &'r Redge<R>) -> Self {
        let handle = redge.vert_handle(vid).edge();
        let hid = if handle.has_hedge() {
            handle.hedge().id()
        } else {
            HedgeId::ABSENT
        };
        Self {
            focused_vertex: vid,
            edge_iter: VertexStarEdgesIter::new(vid, redge),
            orbit_iter: RadialHedgeIter::new(hid, redge),
            seen: BTreeSet::new(),
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for VertexLinkEdgesIter<'r, R> {
    type Item = EdgeHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        // Exhaust the current hedge orbit for all edges in faces touching the vertex.
        if let Some(h) = self.orbit_iter.find(|h| {
            h.source().id() == self.focused_vertex && self.seen.insert(h.face_next().edge().id())
        }) {
            return Some(h.face_next().edge());
        }

        // If the orbit is exhausted, move to the next edge in the star.
        loop {
            let next = self.edge_iter.next();
            if next.is_none() {
                return None;
            }

            let next = next.unwrap();
            if !next.has_hedge() {
                continue;
            }

            self.orbit_iter = RadialHedgeIter::new(next.hedge().id(), self.edge_iter.redge);
            let candidate = self
                .orbit_iter
                .find(|h| {
                    h.source().id() == self.focused_vertex
                        && self.seen.insert(h.face_next().edge().id())
                })
                .map(|h| h.face_next().edge());

            if candidate.is_some() {
                return candidate;
            }
        }
    }
}

pub struct VertIncidentFacesIterator<'r, R: RedgeContainers> {
    edge_iter: VertexStarEdgesIter<'r, R>,
    current_radial_iter: RadialHedgeIter<'r, R>,
    seen: HashSet<FaceId>,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertIncidentFacesIterator<'r, R> {
    pub fn new(vert_id: VertId, redge: &'r Redge<R>) -> Self {
        Self {
            edge_iter: VertexStarEdgesIter::new(vert_id, redge),
            current_radial_iter: RadialHedgeIter::new(
                redge
                    .vert_handle(vert_id)
                    .star_edges()
                    .find(|e| e.has_hedge())
                    .map_or(HedgeId::ABSENT, |e| e.hedge().id()),
                redge,
            ),
            seen: HashSet::new(),
            redge,
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for VertIncidentFacesIterator<'r, R> {
    type Item = FaceHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        // Try to find the next unseen face incident on this edge.
        if let Some(h) = self
            .current_radial_iter
            .find(|h| self.seen.insert(h.face().id()))
        {
            return Some(h.face());
        }

        // If no such face exists, try each edge until we find an unseen incident face or we run out of edges.
        while let Some(e) = self.edge_iter.next() {
            if !e.has_hedge() {
                continue;
            }
            self.current_radial_iter = RadialHedgeIter::new(e.hedge().id(), self.redge);

            if let Some(f) = self
                .current_radial_iter
                .find(|h| self.seen.insert(h.face().id()))
                .map(|h| h.face())
            {
                return Some(f);
            }
        }

        None
    }
}

pub struct VertManifoldIncidentFacesIterator<'r, R: RedgeContainers> {
    last_face: FaceId,
    edge_iter: VertexStarEdgesIter<'r, R>,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertManifoldIncidentFacesIterator<'r, R> {
    pub fn new(vert_id: VertId, redge: &'r Redge<R>) -> Self {
        let vhandle = redge.vert_handle(vert_id);

        Self {
            last_face: FaceId::ABSENT,
            edge_iter: vhandle.star_edges(),
            redge,
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for VertManifoldIncidentFacesIterator<'r, R> {
    type Item = FaceHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let edge = self.edge_iter.next();
            if edge.is_none() {
                return None;
            }
            let edge = edge.unwrap();

            if let Some(fid) = edge
                .hedge()
                .radial_neighbours()
                .map(|h| h.face().id())
                .find(|f| *f != self.last_face)
            {
                self.last_face = fid;
                break;
            }
        }

        Some(FaceHandle::new(self.last_face, self.redge))
    }
}

pub struct RadialHedgeIter<'r, R: RedgeContainers> {
    start_hedge: HedgeId,
    current_hedge: HedgeId,
    start: bool,

    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> RadialHedgeIter<'r, R> {
    pub(crate) fn new(hedge: HedgeId, redge: &'r Redge<R>) -> Self {
        Self {
            start_hedge: hedge,
            current_hedge: hedge,
            start: true,
            redge,
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for RadialHedgeIter<'r, R> {
    type Item = HedgeHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_hedge == HedgeId::ABSENT {
            return None;
        };

        if (self.current_hedge == self.start_hedge && !self.start)
            || self.start_hedge == HedgeId::ABSENT
        {
            return None;
        }

        self.start = false;
        let hedge_handle = self.redge.hedge_handle(self.current_hedge);
        let id = self.current_hedge;
        self.current_hedge = hedge_handle.radial_next().id();

        Some(HedgeHandle::new(id, self.redge))
    }
}

pub struct FaceLoopHedgeIter<'r, R: RedgeContainers> {
    start_hedge: HedgeId,
    current_hedge: HedgeId,
    start: bool,

    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> FaceLoopHedgeIter<'r, R> {
    pub(crate) fn new(hedge: HedgeId, redge: &'r Redge<R>) -> Self {
        Self {
            start_hedge: hedge,
            current_hedge: hedge,
            start: true,
            redge,
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for FaceLoopHedgeIter<'r, R> {
    type Item = HedgeHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_hedge == self.start_hedge && !self.start {
            return None;
        }

        self.start = false;
        let hedge_handle = self.redge.hedge_handle(self.current_hedge);
        let id = self.current_hedge;
        self.current_hedge = hedge_handle.face_next().id();

        Some(HedgeHandle::new(id, self.redge))
    }
}

pub struct FaceVertIterator<'r, R: RedgeContainers> {
    pub(crate) face_loop: FaceLoopHedgeIter<'r, R>,
    pub(crate) _redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> Iterator for FaceVertIterator<'r, R> {
    type Item = VertHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        self.face_loop.next().map(|x| x.source())
    }
}
