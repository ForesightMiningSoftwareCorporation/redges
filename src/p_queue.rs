use std::hash::Hash;

use ordered_float::OrderedFloat;
use priority_queue::PriorityQueue;

/// Max priority queue that can poperly handle floats.
#[derive(Debug)]
pub struct PQueue<T, W>
where
    T: Eq + Hash,
    W: ordered_float::FloatCore,
{
    queue: PriorityQueue<T, OrderedFloat<W>>,
}

impl<T, W> PQueue<T, W>
where
    T: Eq + Hash,
    W: ordered_float::FloatCore,
{
    pub fn new() -> Self {
        PQueue {
            queue: PriorityQueue::new(),
        }
    }

    pub fn push(&mut self, item: T, weight: W) -> Option<OrderedFloat<W>> {
        self.queue.push(item, OrderedFloat(-weight))
    }

    pub fn pop(&mut self) -> Option<(T, W)> {
        self.queue.pop().map(|item| (item.0, item.1 .0))
    }

    pub fn remove(&mut self, item: T) -> Option<(T, W)> {
        self.queue.remove(&item).map(|(i, OrderedFloat(w))| (i, w))
    }

    pub fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }
}
