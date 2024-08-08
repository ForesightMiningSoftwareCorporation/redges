use std::marker::PhantomData;

use crate::{
    container_trait::PrimitiveContainer, container_trait::RedgeContainers, edge_handle::EdgeHandle,
    hedge_handle::HedgeHandle, vert_handle::VertHandle, EdgeId, HedgeId, Redge, VertId,
};

pub struct VertexStarVerticesIter<'r, R: RedgeContainers> {
    start_edge: EdgeId,
    current_edge: EdgeId,
    focused_vertex: VertId,

    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertexStarVerticesIter<'r, R> {
    pub(crate) fn new(vert_id: VertId, redge: &'r Redge<R>) -> Self {
        VertexStarVerticesIter {
            focused_vertex: vert_id,
            start_edge: redge.verts_meta[vert_id.to_index()].edge_id,
            current_edge: redge.verts_meta[vert_id.to_index()].edge_id,
            redge: redge,
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for VertexStarVerticesIter<'r, R> {
    type Item = VertHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        let edge = self.redge.edge_handle(self.current_edge);
        let edge_cycle = edge.edge_cycle_at(self.focused_vertex);

        if edge_cycle.next_edge == self.start_edge {
            return None;
        }

        self.current_edge = edge_cycle.next_edge;

        Some(VertHandle::new(
            edge.opposite(self.focused_vertex),
            self.redge,
        ))
    }
}

pub struct VertexStarEdgesIter<'r, R: RedgeContainers> {
    start_edge: EdgeId,
    current_edge: EdgeId,
    focused_vertex: VertId,

    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertexStarEdgesIter<'r, R> {
    pub(crate) fn new(vert_id: VertId, redge: &'r Redge<R>) -> Self {
        VertexStarEdgesIter {
            focused_vertex: vert_id,
            start_edge: redge.verts_meta[vert_id.to_index()].edge_id,
            current_edge: redge.verts_meta[vert_id.to_index()].edge_id,
            redge: redge,
        }
    }
}

impl<'r, R: RedgeContainers> Iterator for VertexStarEdgesIter<'r, R> {
    type Item = EdgeHandle<'r, R>;

    fn next(&mut self) -> Option<Self::Item> {
        let edge = self.redge.edge_handle(self.current_edge);
        let edge_cycle = edge.edge_cycle_at(self.focused_vertex);

        if edge_cycle.next_edge == self.start_edge {
            return None;
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
        let edge = self.edge_iter.next();

        match edge {
            Some(edge_handle) => Some(EdgeHandle::new(
                edge_handle.hedge().face_next().edge().id(),
                self.edge_iter.redge,
            )),
            None => None,
        }
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
