use crate::{
    container_trait::RedgeContainers, edge_handle::EdgeHandle, hedge_handle::HedgeHandle,
    vert_handle::VertHandle, EdgeId, HedgeId, Redge, VertId,
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
    start_edge: EdgeId,
    current_edge: EdgeId,
    focused_vertex: VertId,
    start: bool,

    redge: &'r Redge<R>,
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

pub struct VertexLinkEdgesIter<'r, R: RedgeContainers> {
    edge_iter: VertexStarEdgesIter<'r, R>,
}

impl<'r, R: RedgeContainers> Iterator for VertexLinkEdgesIter<'r, R> {
    type Item = EdgeHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        self.edge_iter.next().map(|edge_handle| {
            EdgeHandle::new(
                edge_handle.hedge().face_next().edge().id(),
                self.edge_iter.redge,
            )
        })
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
        debug_assert!(self.current_hedge != HedgeId::ABSENT);
        if self.current_hedge == self.start_hedge && !self.start {
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
