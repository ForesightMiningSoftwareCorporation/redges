//! Priority queue that allows for floating point weights.
use std::hash::Hash;

use ordered_float::OrderedFloat;
use priority_queue::PriorityQueue;

/// Max priority queue that can properly handle floats.
#[derive(Debug)]
pub struct PQueue<T, W>
where
    T: Eq + Hash,
    W: ordered_float::FloatCore,
{
    queue: PriorityQueue<T, OrderedFloat<W>>,
}

impl<T, W> Default for PQueue<T, W>
where
    T: Eq + Hash,
    W: ordered_float::FloatCore,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, W> PQueue<T, W>
where
    T: Eq + Hash,
    W: ordered_float::FloatCore,
{
    /// Make a new queue.
    pub fn new() -> Self {
        PQueue {
            queue: PriorityQueue::new(),
        }
    }

    /// Make a queue with pre-allocated size.
    pub fn with_capacity(capacity: usize) -> Self {
        PQueue {
            queue: PriorityQueue::with_capacity(capacity),
        }
    }

    /// Make a queue from an iterator.
    pub fn from_iterator<I>(iter: I) -> Self
    where
        I: Iterator<Item = (T, W)>,
    {
        let mut res = PQueue {
            queue: PriorityQueue::new(),
        };

        for (item, weight) in iter {
            res.push(item, weight);
        }

        res
    }

    /// Insert the `item`, `weight` pair into the queue.
    /// If an element equal to item is already in the queue, its priority is updated and the old
    /// priority is returned in Some; otherwise, item is inserted with priority and `None` is returned.
    pub fn push(&mut self, item: T, weight: W) -> Option<OrderedFloat<W>> {
        self.queue.push(item, OrderedFloat(-weight))
    }

    /// Removes the item with the greatest priority from the priority queue and returns the pair
    /// `item`, `weight`, or `None` if the queue is empty.
    pub fn pop(&mut self) -> Option<(T, W)> {
        if let Some(item) = self.queue.pop() {
            Some((item.0, -item.1 .0))
        } else {
            None
        }
    }

    /// Remove an arbitrary element from the priority queue.
    /// Returns the `item`, `weight` couple or None if the item is not found in the queue.
    /// The operation is performed in O(log(N)) time (worst case).
    pub fn remove(&mut self, item: T) -> Option<(T, W)> {
        self.queue.remove(&item).map(|(i, OrderedFloat(w))| (i, w))
    }
    /// True if the queue has no elements.
    pub fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }
    /// Current number of elements in the queue.
    pub fn len(&self) -> usize {
        self.queue.len()
    }
    /// Iterator over the elements in the queue.
    pub fn iter(&self) -> impl Iterator<Item = (T, W)> + '_
    where
        T: Clone,
    {
        self.queue
            .iter()
            .map(|(item, ord_weight)| (item.clone(), ord_weight.0))
    }
}
