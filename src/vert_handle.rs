use crate::container_trait::{RedgeContainers, VertData};
use crate::edge_handle::EdgeHandle;

use crate::iterators::{
    VertIncidentFacesIterator, VertManifoldIncidentFacesIterator, VertexLinkEdgesIter,
    VertexStarEdgesIter, VertexStarVerticesIter,
};
use crate::EdgeId;
use crate::{container_trait::PrimitiveContainer, Redge, VertId, VertMetaData};

pub struct VertHandle<'r, R: RedgeContainers> {
    id: VertId,
    redge: &'r Redge<R>,
}

impl<'r, R: RedgeContainers> VertHandle<'r, R> {
    pub(crate) fn new(id: VertId, redge: &'r Redge<R>) -> Self {
        debug_assert!(!id.is_absent());
        Self { id, redge }
    }

    pub fn id(&self) -> VertId {
        self.id
    }

    pub fn edge(&self) -> EdgeHandle<'r, R> {
        EdgeHandle::new(
            self.redge.verts_meta[self.id.to_index()].edge_id,
            self.redge,
        )
    }

    pub fn data(&self) -> &VertData<R> {
        self.redge.vert_data.get(self.id.to_index() as u64)
    }

    pub fn is_active(&self) -> bool {
        self.metadata().is_active
    }

    fn metadata(&self) -> &VertMetaData {
        &self.redge.verts_meta[self.id.to_index()]
    }

    pub fn neighbours(&self) -> VertexStarVerticesIter<'r, R> {
        VertexStarVerticesIter::new(self.id, self.redge)
    }

    pub fn star_edges(&self) -> VertexStarEdgesIter<'r, R> {
        VertexStarEdgesIter::new(self.id, self.redge)
    }

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

    pub fn link_edges(&self) -> VertexLinkEdgesIter<'r, R> {
        VertexLinkEdgesIter::new(self.id, self.redge)
    }

    pub fn incident_faces(&self) -> VertIncidentFacesIterator<'r, R> {
        VertIncidentFacesIterator::new(self.id, self.redge)
    }

    /// WARNING: Call only if you know that the vertex is locally manifold.
    /// Otherwise there's no guarantees as to what this method will do.
    pub fn incident_faces_manifold(&self) -> VertManifoldIncidentFacesIterator<'r, R> {
        VertManifoldIncidentFacesIterator::new(self.id, self.redge)
    }

    pub fn is_in_boundary(&self) -> bool {
        self.star_edges().any(|e| e.is_boundary())
    }

    pub fn pick_different(&self, eid: EdgeId) -> Option<EdgeHandle<'r, R>> {
        self.star_edges().find(|e| e.id() != eid)
    }

    pub fn count_incident_boundary_edges(&self) -> usize {
        self.star_edges().filter(|e| e.is_boundary()).count()
    }
}
