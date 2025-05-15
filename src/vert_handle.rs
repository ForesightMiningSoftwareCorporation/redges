//! Topology handle around a vertex.
use crate::container_trait::{RedgeContainers, VertData};
use crate::edge_handle::EdgeHandle;

use crate::iterators::{
    VertIncidentFacesIterator, VertManifoldIncidentFacesIterator, VertexLinkEdgesIter,
    VertexStarEdgesIter, VertexStarVerticesIter,
};
use crate::EdgeId;
use crate::{container_trait::PrimitiveContainer, Redge, VertId, VertMetaData};
/// Topology handle around a vertex.
pub struct VertHandle<'r, R: RedgeContainers> {
    id: VertId,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertHandle<'r, R> {
    pub(crate) fn new(id: VertId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    /// Id.
    pub fn id(&self) -> VertId {
        self.id
    }

    /// Does this vertex point to an edge?
    pub fn has_edge(&self) -> bool {
        self.redge.verts_meta[self.id.to_index()].edge_id != EdgeId::ABSENT
    }
    /// Get the edge. WARNING: Check it exists with `has_edge` if necessary, bad things can happen
    /// otherwise.
    pub fn edge(&self) -> EdgeHandle<'r, R> {
        EdgeHandle::new(
            self.redge.verts_meta[self.id.to_index()].edge_id,
            self.redge,
        )
    }
    /// Get the underlying vertex data.
    pub fn data(&self) -> &VertData<R> {
        self.redge.vert_data.get(self.id.to_index() as u64)
    }
    /// Is this vertex active/valid?
    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &VertMetaData {
        &self.redge.verts_meta[self.id.to_index()]
    }
    /// Iterator around all vertices sharing an edge with this vertex.
    pub fn neighbours(&self) -> VertexStarVerticesIter<'r, R> {
        VertexStarVerticesIter::new(self.id, self.redge)
    }
    /// Iterator around all edges incident on the vertex *in no particular order*. Do *NOT* assume
    /// this has a winding order.
    pub fn star_edges(&self) -> VertexStarEdgesIter<'r, R> {
        VertexStarEdgesIter::new(self.id, self.redge)
    }
    /// Like `star_edges` but starts the iteration at the specified edge.
    /// Assumes the edge is in the star.
    pub fn star_edges_from(&self, edge_id: EdgeId) -> VertexStarEdgesIter<'r, R> {
        debug_assert!(
            self.redge.edge_handle(edge_id).v1().id() == self.id
                || self.redge.edge_handle(edge_id).v2().id() == self.id
        );
        VertexStarEdgesIter {
            start_edge: edge_id,
            current_edge: edge_id,
            focused_vertex: self.id,
            start: true,
            redge: self.redge,
        }
    }
    /// Iterator around the edges of the link of the vertex.
    pub fn link_edges(&self) -> VertexLinkEdgesIter<'r, R> {
        VertexLinkEdgesIter::new(self.id, self.redge)
    }
    /// Iterator over the faces incident on the vertex.
    pub fn incident_faces(&self) -> VertIncidentFacesIterator<'r, R> {
        VertIncidentFacesIterator::new(self.id, self.redge)
    }

    /// WARNING: Call only if you know that the vertex is locally manifold.
    /// Otherwise, there are no guarantees as to what this method will do.
    pub fn incident_faces_manifold(&self) -> VertManifoldIncidentFacesIterator<'r, R> {
        VertManifoldIncidentFacesIterator::new(self.id, self.redge)
    }

    /// Is this vertex pointed to by a boundary edge?
    pub fn is_in_boundary(&self) -> bool {
        self.star_edges().any(|e| e.is_boundary())
    }

    /// Pick an edge in the star other than the specified one, if possible.
    pub fn pick_different(&self, eid: EdgeId) -> Option<EdgeHandle<'r, R>> {
        self.star_edges().find(|e| e.id() != eid)
    }

    /// Count the number of edges in the boundary that share this vertex.
    pub fn count_incident_boundary_edges(&self) -> usize {
        self.star_edges().filter(|e| e.is_boundary()).count()
    }
}
