use std::{fmt::Debug, usize};

use linear_isomorphic::{InnerSpace, RealField};
use nalgebra::Vector3;

use crate::VertId;

pub trait PrimitiveContainer: Clone + Debug {
    type PrimitiveData: Clone + Debug;

    fn get(&self, index: u64) -> &Self::PrimitiveData;
    fn get_mut(&mut self, index: u64) -> &mut Self::PrimitiveData;
    fn set(&mut self, index: u64, data: Self::PrimitiveData);

    fn push(&mut self, data: Self::PrimitiveData);
    fn remove(&mut self, index: u64) -> Self::PrimitiveData;
    fn swap_remove(&mut self, index: u64) -> Self::PrimitiveData;

    fn retain<F: FnMut(&Self::PrimitiveData) -> bool>(&mut self, f: F);

    fn resize(&mut self, new_size: usize);

    fn len(&self) -> usize;

    fn iterate(&self) -> impl Iterator<Item = &Self::PrimitiveData>;
}

impl<T: Default + Clone + Debug> PrimitiveContainer for Vec<T> {
    type PrimitiveData = T;

    fn get(&self, index: u64) -> &Self::PrimitiveData {
        &self[index as usize]
    }

    fn get_mut(&mut self, index: u64) -> &mut Self::PrimitiveData {
        &mut self[index as usize]
    }

    fn set(&mut self, index: u64, data: Self::PrimitiveData) {
        self[index as usize] = data
    }

    fn push(&mut self, data: Self::PrimitiveData) {
        self.push(data)
    }

    fn remove(&mut self, index: u64) -> Self::PrimitiveData {
        self.remove(index as usize)
    }

    fn swap_remove(&mut self, index: u64) -> Self::PrimitiveData {
        self.swap_remove(index as usize)
    }

    fn retain<F: FnMut(&Self::PrimitiveData) -> bool>(&mut self, f: F) {
        self.retain(f)
    }

    fn resize(&mut self, new_size: usize) {
        self.resize(new_size, T::default());
    }

    fn iterate(&self) -> impl Iterator<Item = &Self::PrimitiveData> {
        self.iter()
    }

    fn len(&self) -> usize {
        self.len()
    }
}

impl PrimitiveContainer for () {
    type PrimitiveData = ();

    fn get(&self, _index: u64) -> &Self::PrimitiveData {
        self
    }

    fn get_mut(&mut self, _index: u64) -> &mut Self::PrimitiveData {
        self
    }

    fn set(&mut self, _index: u64, _data: Self::PrimitiveData) {}

    fn push(&mut self, _he1data: Self::PrimitiveData) {}

    fn remove(&mut self, _index: u64) -> Self::PrimitiveData {}

    fn swap_remove(&mut self, _index: u64) -> Self::PrimitiveData {
        ()
    }

    fn resize(&mut self, _new_size: usize) {}

    fn retain<F: FnMut(&Self::PrimitiveData) -> bool>(&mut self, _f: F) {}

    fn iterate(&self) -> impl Iterator<Item = &Self::PrimitiveData> {
        std::iter::empty()
    }

    fn len(&self) -> usize {
        0
    }
}

pub trait RedgeContainers {
    type VertContainer: PrimitiveContainer;
    type EdgeContainer: PrimitiveContainer;
    type FaceContainer: PrimitiveContainer;
}

pub type VertData<R> = <<R as RedgeContainers>::VertContainer as PrimitiveContainer>::PrimitiveData;
pub type EdgeData<R> = <<R as RedgeContainers>::EdgeContainer as PrimitiveContainer>::PrimitiveData;
pub type FaceData<R> = <<R as RedgeContainers>::FaceContainer as PrimitiveContainer>::PrimitiveData;

impl<V, E, F> RedgeContainers for (V, E, F)
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
    type VertContainer = V;
    type EdgeContainer = E;
    type FaceContainer = F;
}

pub trait FaceAttributeGetter<S> {
    fn attribute_count(&self) -> usize;
    fn attribute(&self, vert_index: usize, attribute_id: usize) -> S;
    fn attribute_mut(&mut self, vert_index: usize, attribute_id: usize) -> &mut S;
    fn inner_index(&self, vid: VertId) -> usize;
    fn attribute_vertices(&self) -> &[VertId];
    fn attribute_vertices_mut(&mut self) -> &mut [VertId];
}

impl<S> FaceAttributeGetter<S> for () {
    fn attribute_count(&self) -> usize {
        0
    }

    fn attribute(&self, _vert_index: usize, _attribute_id: usize) -> S {
        panic!()
    }

    fn attribute_mut(&mut self, _vert_index: usize, _attribute_id: usize) -> &mut S {
        panic!()
    }

    fn attribute_vertices(&self) -> &[VertId] {
        panic!()
    }

    fn attribute_vertices_mut(&mut self) -> &mut [VertId] {
        panic!()
    }

    fn inner_index(&self, _vid: VertId) -> usize {
        usize::MAX
    }
}

/// The `x,y,z` coordinates should be attributes 0, 1, 2 respectively. So
/// any vertex *must* have at least three attributes.
pub trait VertexAttributeGetter<S> {
    fn attribute_count(&self) -> usize;
    fn attribute(&self, attribute_id: usize) -> S;
}

impl<S> VertexAttributeGetter<S> for Vector3<S>
where
    S: RealField,
{
    fn attribute_count(&self) -> usize {
        3
    }

    fn attribute(&self, attribute_id: usize) -> S {
        self[attribute_id]
    }
}
