// The functions defined here should remain functions, don't refactor them to method.
// The reason is that, we will likely want to implement trait interfaces
// to abstract meshes. These methods should work for those abstractions.
// If we couple them to the objects it will be harder to do this.

use crate::{PrimitiveContainer, Redge};

pub fn is_manifold<V, E, F>(mesh: &Redge<V, E, F>)
where
    V: PrimitiveContainer,
    E: PrimitiveContainer,
    F: PrimitiveContainer,
{
}
