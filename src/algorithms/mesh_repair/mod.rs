//! Mesh repair algorithms for fixing common defects.
//!
//! This module provides tools for repairing mesh topology and geometry:
//! - Vertex welding: Merge duplicate or nearby vertices
//! - Future: Non-manifold repair, hole filling, degeneracy removal

pub mod vertex_welding;

pub use vertex_welding::VertexWeldingReport;
