use std::marker::PhantomData;

use crate::{edge_handle, vert_handle::VertHandle, EdgeId, PrimitiveContainer, Redge, VertId};

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
    type Item = VertId;

    fn next(&mut self) -> Option<Self::Item> {
        let edge = self.redge.edge_handle(self.current_edge);
        let edge_cycle = edge.edge_cycle_at(self.focused_vertex);

        if edge_cycle.next_edge == self.start_edge {
            return None;
        }

        self.current_edge = edge_cycle.next_edge;

        Some(edge.opposite(self.focused_vertex))
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
    type Item = EdgeId;

    fn next(&mut self) -> Option<Self::Item> {
        let edge = self.redge.edge_handle(self.current_edge);
        let edge_cycle = edge.edge_cycle_at(self.focused_vertex);

        if edge_cycle.next_edge == self.start_edge {
            return None;
        }

        self.current_edge = edge_cycle.next_edge;

        Some(edge.id())
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
    type Item = EdgeId;

    fn next(&mut self) -> Option<Self::Item> {
        let edge = self.edge_iter.next();

        match edge {
            Some(edge_id) => {
                let edge_handle = self.edge_iter.redge.edge_handle(edge_id);
                Some(edge_handle.hedge().face_next().edge().id())
            }
            None => None,
        }
    }
}
