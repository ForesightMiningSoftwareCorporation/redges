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

    pub fn with_capacity(capacity: usize) -> Self {
        PQueue {
            queue: PriorityQueue::with_capacity(capacity),
        }
    }

    pub fn from_iter<I>(iter: I) -> Self
    where
        I: Iterator<Item = (T, W)>,
    {
        let mut res = PQueue {
            queue: PriorityQueue::new(),
        };

        for (item, weight) in iter {
            res.push(-weight, item);
        }

        res
    }

    pub fn push(&mut self, weight: W, item: T) -> Option<OrderedFloat<W>> {
        self.queue.push(item, OrderedFloat(-weight))
    }

    pub fn pop(&mut self) -> Option<(W, T)> {
        if let Some(item) = self.queue.pop() {
            Some((-item.1 .0, item.0))
        } else {
            None
        }
    }

    pub fn remove(&mut self, item: T) -> Option<(T, W)> {
        match self.queue.remove(&item) {
            None => None,
            Some((i, OrderedFloat(w))) => Some((i, w)),
        }
    }

    pub fn contains(&self, item: &T) -> bool {
        self.queue.get(item).is_some()
    }

    pub fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }

    pub fn len(&self) -> usize {
        self.queue.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = (T, W)> + '_
    where
        T: Clone,
    {
        self.queue
            .iter()
            .map(|(item, ord_weight)| (item.clone(), ord_weight.0))
    }
}
