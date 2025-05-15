//! Traits to abstract over a collection/container.

use std::fmt::Debug;

use linear_isomorphic::RealField;
use nalgebra::Vector3;

use crate::VertId;

/// Trait to abstract over a collection, usually a `Vec<T>`.
#[allow(clippy::len_without_is_empty)]
pub trait PrimitiveContainer: Clone + Debug {
    /// Underlying datum of the container.
    type PrimitiveData: Clone + Debug;

    /// Get datum.
    fn get(&self, index: u64) -> &Self::PrimitiveData;
    /// Get mutable datum.
    fn get_mut(&mut self, index: u64) -> &mut Self::PrimitiveData;
    /// Set datum.
    fn set(&mut self, index: u64, data: Self::PrimitiveData);

    /// Add datum to container.
    fn push(&mut self, data: Self::PrimitiveData);
    /// Remove element at `index` from container.
    fn remove(&mut self, index: u64) -> Self::PrimitiveData;
    /// Remove element at `index` from container and replace it with the last element.
    fn swap_remove(&mut self, index: u64) -> Self::PrimitiveData;
    /// Retain a subset of the container, delete the rest.
    fn retain<F: FnMut(&Self::PrimitiveData) -> bool>(&mut self, f: F);
    /// Increment or reduce the size of the container.
    fn resize(&mut self, new_size: usize);
    /// Cardinality of the container.
    fn len(&self) -> usize;
    /// Iterate over the container data.
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

    fn len(&self) -> usize {
        self.len()
    }

    fn iterate(&self) -> impl Iterator<Item = &Self::PrimitiveData> {
        self.iter()
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

    fn swap_remove(&mut self, _index: u64) -> Self::PrimitiveData {}

    fn resize(&mut self, _new_size: usize) {}

    fn retain<F: FnMut(&Self::PrimitiveData) -> bool>(&mut self, _f: F) {}

    fn iterate(&self) -> impl Iterator<Item = &Self::PrimitiveData> {
        std::iter::empty()
    }

    fn len(&self) -> usize {
        0
    }
}

/// Helper trait to define the geometry information fo a radial edge.
pub trait RedgeContainers {
    /// Vertex data container.
    type VertContainer: PrimitiveContainer;
    /// Edge data container.
    type EdgeContainer: PrimitiveContainer;
    /// Vertex data container.
    type FaceContainer: PrimitiveContainer;
}

/// Underlying data of a vertex.
pub type VertData<R> = <<R as RedgeContainers>::VertContainer as PrimitiveContainer>::PrimitiveData;
/// Underlying data of an edge.
pub type EdgeData<R> = <<R as RedgeContainers>::EdgeContainer as PrimitiveContainer>::PrimitiveData;
/// Underlying data of a face.
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

/// Trait to enable querying information about face data, such as normals or uvs.
pub trait FaceAttributeGetter<S> {
    /// Number of *scalar* attributes. So for example, a normal
    /// adds *three* attributes and a uv adds *two*. So
    /// a face with both normals and uvs has five attributes.
    fn attribute_count(&self) -> usize;
    /// Fetch a *scalar* attribute, the trait assumes consistent ordering.
    fn attribute(&self, vert_index: usize, attribute_id: usize) -> S;
    /// Fetch a *scalar* attribute, the trait assumes consistent ordering.
    fn attribute_mut(&mut self, vert_index: usize, attribute_id: usize) -> &mut S;
    /// Return the local index of a vertex in this face, e.g. in a triangle mesh
    /// the possible values are 0, 1, 2.
    fn inner_index(&self, vid: VertId) -> usize;
    /// Constant ordered list of vertex ids in this face.
    /// This list must maintain the same order, regardless of how the
    /// mesh is structured.
    // TODO: Might be enough to stort them by id?
    fn attribute_vertices(&self) -> &[VertId];
    /// Similar as `attribute_vertices` but allows for updating.
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
    /// Number of scalar attribute values in this vertex.
    /// This is at least the number of coordinates in the vertex, normally three.
    fn attribute_count(&self) -> usize;
    /// Scalar attribute value at the given attribute index.
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
