//! Unified mesh repair orchestrator and configuration.
//!
//! This module provides a high-level API for mesh repair that orchestrates
//! multiple repair algorithms in the correct order. Individual algorithms can
//! also be used standalone via their respective modules.
//!
//! # Architecture
//!
//! - **Individual algorithms**: Each in its own module (vertex_welding, non_manifold_repair, etc.)
//! - **This module**: Orchestrates algorithms via `MeshRepairConfig` and `repair_mesh()`
//! - **Re-exports**: All individual algorithms are re-exported for standalone use
//!
//! # Pipeline Order
//!
//! Repair operations must be performed in a specific order because some operations
//! can create or expose new defects:
//!
//! 1. **Vertex welding** - Merge duplicate vertices (can create degenerate faces)
//! 2. **Remove degeneracies** - Clean up zero-area faces/edges
//! 3. **Non-manifold repair** - Split non-manifold edges and vertices
//! 4. **Hole filling** - Fill boundary loops (future)
//! 5. **Normal orientation** - Ensure consistent face orientation (future)
//!
//! The pipeline can iterate multiple times if early operations expose new defects.
//!
//! # Example: Full Repair Pipeline
//!
//! ```rust,ignore
//! use redges::algorithms::mesh_repair::{MeshRepairConfig, repair_mesh};
//!
//! let mesh = load_mesh("broken_model.obj");
//!
//! let config = MeshRepairConfig {
//!     weld_tolerance: Some(1e-6),
//!     remove_degeneracies: true,
//!     fix_nonmanifold: true,
//!     max_iterations: 3,
//!     ..Default::default()
//! };
//!
//! let (repaired_mesh, report) = mesh.repair_mesh(config);
//!
//! println!("Repair complete in {} iterations", report.iterations_needed);
//! println!("  Vertices merged: {}", report.vertices_merged);
//! println!("  Degenerate faces removed: {}", report.degenerate_faces_removed);
//! println!("  Non-manifold elements fixed: {}", report.nonmanifold_elements_fixed);
//! ```
//!
//! # Example: Individual Algorithm
//!
//! ```rust,ignore
//! use redges::algorithms::mesh_repair::vertex_welding::weld_vertices;
//!
//! let mesh = load_mesh("model.obj");
//! let (welded, report) = mesh.weld_vertices(1e-6);
//! ```

// Re-export individual repair algorithms for standalone use
pub mod non_manifold_repair;
pub mod vertex_welding;

// Re-export commonly used types from individual modules
pub use non_manifold_repair::NonManifoldRepairReport;
pub use vertex_welding::VertexWeldingReport;

use crate::{container_trait::RedgeContainers, mesh_deleter::MeshDeleter, Redge};
use linear_isomorphic::{InnerSpace, RealField};
use num::{Bounded, Signed};

/// Configuration for unified mesh repair pipeline.
///
/// This allows fine-grained control over which repair operations are performed
/// and in what order. Operations are always performed in the optimal sequence
/// regardless of configuration order.
#[derive(Debug, Clone)]
pub struct MeshRepairConfig<S> {
    /// Tolerance for vertex welding. If Some(tolerance), weld vertices within this distance.
    /// If None, skip vertex welding.
    ///
    /// Recommended values:
    /// - 3D printing: 0.01mm to 0.1mm
    /// - CAD models: 1e-6 to 1e-3
    /// - Exact duplicates only: 1e-10
    pub weld_tolerance: Option<S>,

    /// Whether to remove degenerate faces and edges.
    ///
    /// Degenerate geometry includes:
    /// - Zero-area faces (collinear vertices)
    /// - Zero-length edges (coincident endpoints)
    /// - Faces with duplicate vertices
    ///
    /// This is typically enabled because degeneracies cause problems
    /// for most downstream operations.
    pub remove_degeneracies: bool,

    /// Whether to repair non-manifold topology.
    ///
    /// This splits:
    /// - Non-manifold edges (edges with 3+ faces)
    /// - Non-manifold vertices (bowtie configurations)
    ///
    /// Note: Currently only detects non-manifold elements (repair is placeholder).
    pub fix_nonmanifold: bool,

    /// Whether to fill holes (boundary loops).
    ///
    /// This creates new faces to close gaps in the mesh, making it watertight.
    ///
    /// Note: Not yet implemented.
    pub fill_holes: bool,

    /// Whether to orient face normals consistently.
    ///
    /// This ensures all faces point "outward" from the solid interior.
    ///
    /// Note: Not yet implemented.
    pub orient_normals: bool,

    /// Whether to detect and repair self-intersections.
    ///
    /// This is expensive and typically opt-in.
    ///
    /// Note: Not yet implemented.
    pub fix_intersections: bool,

    /// Maximum number of repair iterations.
    ///
    /// The pipeline can iterate multiple times because:
    /// - Vertex welding can create degenerate faces
    /// - Degeneracy removal can create non-manifold topology
    /// - Hole filling can create new defects
    ///
    /// Typical values: 1-3 iterations
    pub max_iterations: usize,
}

impl<S> Default for MeshRepairConfig<S>
where
    S: From<f64>,
{
    fn default() -> Self {
        Self {
            weld_tolerance: Some(S::from(1e-6)),
            remove_degeneracies: true,
            fix_nonmanifold: true,
            fill_holes: false,         // Not implemented yet
            orient_normals: false,      // Not implemented yet
            fix_intersections: false,   // Not implemented yet, expensive
            max_iterations: 3,
        }
    }
}

/// Report from unified mesh repair pipeline.
///
/// This aggregates reports from all individual repair operations,
/// providing a comprehensive view of what was fixed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MeshRepairReport {
    /// Number of vertices that were merged during vertex welding.
    pub vertices_merged: usize,

    /// Number of degenerate faces removed.
    pub degenerate_faces_removed: usize,

    /// Number of degenerate edges removed.
    pub degenerate_edges_removed: usize,

    /// Number of non-manifold edges split.
    pub edges_split: usize,

    /// Number of non-manifold vertices split.
    pub vertices_split: usize,

    /// Number of holes filled (not yet implemented).
    pub holes_filled: usize,

    /// Number of faces flipped for consistent orientation (not yet implemented).
    pub faces_flipped: usize,

    /// Number of self-intersections repaired (not yet implemented).
    pub intersections_fixed: usize,

    /// Number of iterations needed to stabilize the mesh.
    pub iterations_needed: usize,
}

impl Default for MeshRepairReport {
    fn default() -> Self {
        Self {
            vertices_merged: 0,
            degenerate_faces_removed: 0,
            degenerate_edges_removed: 0,
            edges_split: 0,
            vertices_split: 0,
            holes_filled: 0,
            faces_flipped: 0,
            intersections_fixed: 0,
            iterations_needed: 0,
        }
    }
}

impl MeshRepairReport {
    /// Compute total number of non-manifold elements fixed.
    pub fn nonmanifold_elements_fixed(&self) -> usize {
        self.edges_split + self.vertices_split
    }

    /// Check if any repairs were performed.
    pub fn any_repairs(&self) -> bool {
        self.vertices_merged > 0
            || self.degenerate_faces_removed > 0
            || self.degenerate_edges_removed > 0
            || self.edges_split > 0
            || self.vertices_split > 0
            || self.holes_filled > 0
            || self.faces_flipped > 0
            || self.intersections_fixed > 0
    }
}

impl<R: RedgeContainers> Redge<R> {
    /// Repair mesh using the specified configuration.
    ///
    /// This orchestrates multiple repair operations in the optimal order:
    /// 1. Vertex welding (if configured)
    /// 2. Degeneracy removal (if configured)
    /// 3. Non-manifold repair (if configured)
    /// 4. Hole filling (if configured, not yet implemented)
    /// 5. Normal orientation (if configured, not yet implemented)
    ///
    /// The pipeline can iterate multiple times if operations expose new defects.
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration specifying which operations to perform
    ///
    /// # Returns
    ///
    /// Returns the repaired mesh and a detailed report of changes made.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let config = MeshRepairConfig::default();
    /// let (repaired, report) = mesh.repair_mesh(config);
    /// ```
    pub fn repair_mesh<S>(mut self, config: MeshRepairConfig<S>) -> (Self, MeshRepairReport)
    where
        crate::container_trait::VertData<R>: InnerSpace<S>,
        S: RealField + Bounded + Signed + Copy,
    {
        let mut total_report = MeshRepairReport::default();
        total_report.iterations_needed = 1; // At minimum, we always run one iteration

        for iteration in 0..config.max_iterations {
            let mut iteration_report = MeshRepairReport::default();

            // Step 1: Vertex Welding
            if let Some(tolerance) = config.weld_tolerance {
                let (welded_mesh, weld_report) = self.weld_vertices(tolerance);
                self = welded_mesh;

                iteration_report.vertices_merged += weld_report.vertices_merged;
                iteration_report.degenerate_faces_removed += weld_report.degenerate_faces_removed;
                iteration_report.degenerate_edges_removed += weld_report.degenerate_edges_removed;
            }

            // Step 2: Degeneracy Removal
            if config.remove_degeneracies {
                // TODO: Implement comprehensive degeneracy removal
                // For now, vertex welding handles some degeneracies
            }

            // Step 3: Non-Manifold Repair
            if config.fix_nonmanifold {
                let mut deleter = MeshDeleter::start_deletion(self);
                let nm_report = deleter.repair_nonmanifold_topology();

                // Only count as a repair if something was actually split
                if nm_report.edges_split > 0 || nm_report.vertices_split > 0 {
                    iteration_report.edges_split += nm_report.edges_split;
                    iteration_report.vertices_split += nm_report.vertices_split;
                }

                self = deleter.end_deletion();
            }

            // Step 4: Hole Filling (future)
            if config.fill_holes {
                // TODO: Implement hole filling
            }

            // Step 5: Normal Orientation (future)
            if config.orient_normals {
                // TODO: Implement normal orientation
            }

            // Step 6: Self-Intersection Repair (future)
            if config.fix_intersections {
                // TODO: Implement intersection repair
            }

            // If nothing changed this iteration, we're done
            if !iteration_report.any_repairs() {
                break;
            }

            // Accumulate reports (only if we did something)
            total_report.vertices_merged += iteration_report.vertices_merged;
            total_report.degenerate_faces_removed += iteration_report.degenerate_faces_removed;
            total_report.degenerate_edges_removed += iteration_report.degenerate_edges_removed;
            total_report.edges_split += iteration_report.edges_split;
            total_report.vertices_split += iteration_report.vertices_split;
            total_report.iterations_needed = iteration + 1;
        }

        (self, total_report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_repair_mesh_with_default_config() {
        // Create a mesh with duplicate vertices
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0), // Duplicate
        ];
        let faces = vec![vec![0, 1, 2]];
        let mesh = Redge::<(Vec<Vector3<f64>>, (), ())>::new(
            vertices,
            (),
            (),
            faces.into_iter().map(|f| f.into_iter()),
        );

        let config = MeshRepairConfig::default();
        let (_repaired, report) = mesh.repair_mesh(config);

        // Should have merged the duplicate vertex
        assert!(report.vertices_merged > 0);
        assert_eq!(report.iterations_needed, 1);
    }

    #[test]
    fn test_repair_mesh_skip_vertex_welding() {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0), // Duplicate
        ];
        let faces = vec![vec![0, 1, 2]];
        let mesh = Redge::<(Vec<Vector3<f64>>, (), ())>::new(
            vertices,
            (),
            (),
            faces.into_iter().map(|f| f.into_iter()),
        );

        let config = MeshRepairConfig {
            weld_tolerance: None, // Skip welding
            ..Default::default()
        };
        let (_repaired, report) = mesh.repair_mesh(config);

        // Should not have merged vertices
        assert_eq!(report.vertices_merged, 0);
    }

    #[test]
    fn test_repair_clean_mesh_is_noop() {
        // Clean mesh with no defects
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, 1.0, 0.0),
        ];
        let faces = vec![vec![0, 1, 2]];
        let mesh = Redge::<(Vec<Vector3<f64>>, (), ())>::new(
            vertices,
            (),
            (),
            faces.into_iter().map(|f| f.into_iter()),
        );

        let config = MeshRepairConfig::default();
        let (_repaired, report) = mesh.repair_mesh(config);

        // Should not have changed anything
        assert!(!report.any_repairs());
        assert_eq!(report.iterations_needed, 1);
    }
}
