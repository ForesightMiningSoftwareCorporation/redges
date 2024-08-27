use std::fmt::Debug;

/// An element's position in the `heap`.
pub type HeapIndex = u32;

/// An element's position in the `values` vector.
pub type ValueIndex = u32;

/// Binary heap that tracks the indices of its elements.
///
/// This supports accessing specific elements by an independent index.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct IndexBinaryHeap<T> {
    pub values: Vec<T>,
    pub heap: Vec<u32>,
    pub indices: Vec<HeapIndex>,
}

impl<T: Copy + PartialOrd + Default + Debug> IndexBinaryHeap<T> {
    // create a new instance of the heap
    pub fn new() -> Self {
        IndexBinaryHeap {
            values: Vec::new(),
            heap: Vec::new(),
            indices: Vec::new(),
        }
    }

    pub fn reserve(&mut self, additional: usize) {
        self.heap.reserve(additional);
        self.indices.reserve(additional);
        self.values.reserve(additional);
    }

    pub fn peek(&self) -> Option<(T, u32)> {
        if self.heap.is_empty() {
            return None;
        }

        Some((self.values[self.heap[0] as usize], self.heap[0]))
    }

    // get the length of the heap
    pub fn len(&self) -> usize {
        self.heap.len()
    }

    pub fn clear(&mut self) {
        self.heap.clear();
        self.indices.clear();
        self.values.clear();
    }

    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    // push a new element onto the heap
    pub fn push(&mut self, item: T, index: ValueIndex) {
        let old_len = self.heap.len();
        let new_indices_len = self.indices.len().max(index as usize + 1);
        self.indices.resize(new_indices_len, u32::MAX);
        self.values.resize(new_indices_len, Default::default());

        self.indices[index as usize] = old_len as u32;
        self.values[index as usize] = item;
        self.heap.push(index);

        self.sift_up(old_len as u32);
    }

    // pop the element with the highest value from the heap
    pub fn pop(&mut self) -> Option<(T, ValueIndex)> {
        if self.heap.is_empty() {
            return None;
        }

        let idx = self.heap[0];
        let value = self.values[self.heap[0] as usize];

        self.indices[idx as usize] = u32::MAX;
        self.heap.swap_remove(0);

        if !self.heap.is_empty() {
            self.indices[self.heap[0] as usize] = u32::MAX;
            self.sift_down(0);
        }

        Some((value, idx))
    }

    pub fn contains(&self, index: ValueIndex) -> bool {
        index < self.indices.len() as u32
            && self.indices[index as usize] != u32::MAX
            && self.indices[index as usize] < self.heap.len() as u32
    }

    pub fn remove(&mut self, index: ValueIndex) {
        if !self.contains(index) {
            return;
        }

        let value = self.values[index as usize];
        let heap_idx = self.indices[index as usize];

        self.heap[heap_idx as usize] = *self.heap.last().unwrap();
        self.indices[self.heap[heap_idx as usize] as usize] = heap_idx;
        self.indices[index as usize] = u32::MAX;

        if value < self.values[self.heap[heap_idx as usize] as usize] {
            self.sift_down(heap_idx);
        } else {
            self.sift_up(heap_idx);
        }

        self.heap.pop();
    }

    #[inline]
    pub fn update(&mut self, item: T, index: ValueIndex) {
        self.values[index as usize] = item;
        let heap_idx = self.indices[index as usize];

        if heap_idx > 0 && item < self.values[self.heap[((heap_idx - 1) / 2) as usize] as usize] {
            self.sift_up(heap_idx);
        } else {
            self.sift_down(heap_idx);
        }
    }

    fn sift_up(&mut self, index: HeapIndex) {
        let moved = self.heap[index as usize];
        let mut i = index as usize;

        while i > 0 {
            let parent = (i - 1) / 2;

            if self.values[moved as usize] >= self.values[self.heap[parent] as usize] {
                break;
            }

            self.heap[i] = self.heap[parent];
            self.indices[self.heap[i] as usize] = i as u32;

            i = parent;
        }

        if i != index as usize {
            self.heap[i] = moved;
            self.indices[moved as usize] = i as u32;
        }
    }

    fn sift_down(&mut self, index: HeapIndex) {
        let moved = self.heap[index as usize];
        let mut i = index as usize;
        let mut left = (i * 2) + 1;
        let mut right = left + 1;

        while left < self.heap.len() {
            let mut smallest = left;
            if right < self.heap.len()
                && self.values[self.heap[left] as usize] >= self.values[self.heap[right] as usize]
            {
                smallest = right;
            }

            if self.values[self.heap[smallest] as usize] < self.values[moved as usize] {
                self.heap[i] = self.heap[smallest];
                self.indices[self.heap[i] as usize] = i as u32;

                i = smallest;
                left = (i * 2) + 1;
                right = left + 1;
            } else {
                break;
            }
        }

        if i != index as usize {
            self.heap[i] = moved;
            self.indices[moved as usize] = i as u32;
        }
    }
}

#[test]
fn binary_heap_test() {
    let mut heap = IndexBinaryHeap::new();
    heap.push(1, 0);
    heap.push(2, 1);
    heap.push(3, 2);
    heap.push(4, 3);
    heap.push(5, 4);
    heap.push(6, 5);
    heap.push(7, 6);
    heap.push(8, 7);
    heap.push(9, 8);
    heap.push(10, 9);

    for (i, idx) in heap.indices.iter().enumerate() {
        assert_eq!(i as u32, heap.heap[*idx as usize]);
    }

    assert_eq!(heap.len(), 10);

    heap.update(4, 9);

    heap.remove(4);
    heap.remove(8);

    assert_eq!(heap.pop(), Some((1, 0)));
    assert_eq!(heap.pop(), Some((2, 1)));
    assert_eq!(heap.pop(), Some((3, 2)));
    assert_eq!(heap.pop(), Some((4, 9)));

    assert_eq!(heap.pop(), Some((4, 3)));

    assert_eq!(heap.pop(), Some((6, 5)));
    assert_eq!(heap.pop(), Some((7, 6)));
    assert_eq!(heap.pop(), Some((8, 7)));

    assert_eq!(heap.len(), 0);

    heap.push(1, 0);

    assert_eq!(heap.len(), 1);

    assert_eq!(heap.pop(), Some((1, 0)));
}
