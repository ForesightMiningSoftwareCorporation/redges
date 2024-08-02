use std::marker::PhantomData;

use crate::{
    edge_handle::EdgeHandle, hedge_handle::HedgeHandle, vert_handle::VertHandle, EdgeId, HedgeId,
    PrimitiveContainer, Redge, VertId,
};

pub struct VertexStarVerticesIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    start_edge: EdgeId,
    current_edge: EdgeId,
    focused_vertex: VertId,

    redge: &'r Redge<V, E, F>,
}

impl<'r, V, E, F> Iterator for VertexStarVerticesIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    type Item = VertHandle<'r, V, E, F>;

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

pub struct VertexStarEdgesIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    start_edge: EdgeId,
    current_edge: EdgeId,
    focused_vertex: VertId,

    redge: &'r Redge<V, E, F>,
}

impl<'r, V, E, F> Iterator for VertexStarEdgesIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    type Item = EdgeHandle<'r, V, E, F>;

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

pub struct VertexLinkEdgesIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    edge_iter: VertexStarEdgesIter<'r, V, E, F>,
}

impl<'r, V, E, F> Iterator for VertexLinkEdgesIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    type Item = EdgeHandle<'r, V, E, F>;

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

pub struct RadialHedgeIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    start_hedge: HedgeId,
    current_hedge: HedgeId,
    start: bool,

    redge: &'r Redge<V, E, F>,
}

impl<'r, V, E, F> RadialHedgeIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    pub(crate) fn new(hedge: HedgeId, redge: &'r Redge<V, E, F>) -> Self {
        Self {
            start_hedge: hedge,
            current_hedge: hedge,
            start: true,
            redge,
        }
    }
}

impl<'r, V, E, F> Iterator for RadialHedgeIter<'r, V, E, F>
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    type Item = HedgeHandle<'r, V, E, F>;

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
