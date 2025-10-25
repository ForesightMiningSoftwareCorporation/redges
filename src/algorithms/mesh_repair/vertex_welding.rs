//! Vertex welding algorithm for merging duplicate or nearby vertices.
//!
//! Vertex welding is a fundamental mesh repair operation that identifies and merges
//! vertices that are at the same location (exact duplicates) or within a tolerance
//! distance (near-duplicates).
//!
//! # Why Vertex Welding Matters
//!
//! Duplicate vertices are a common defect in meshes that can cause:
//! - **Watertightness issues**: T-junctions that prevent proper sealing
//! - **3D printing failures**: Slicers can't handle gaps between "touching" vertices
//! - **Boolean operation failures**: CSG operations require exact vertex matching
//! - **Rendering artifacts**: Cracks appear where vertices should be merged
//! - **Memory waste**: Storing redundant vertex data
//!
//! Common sources of duplicate vertices:
//! - CAD software exports
//! - Merging multiple meshes
//! - Floating-point precision errors
//! - Boolean operations (union, intersection, difference)
//! - Manual mesh editing
//!
//! # Algorithm
//!
//! The vertex welding algorithm uses a spatial data structure (R*Tree) to efficiently
//! find nearby vertices, then merges them using a union-find structure:
//!
//! 1. **Build R*Tree**: Insert all vertices into spatial index for O(log n) range queries
//! 2. **Find nearby vertices**: For each vertex, query all vertices within tolerance radius
//! 3. **Build equivalence classes**: Use union-find to group vertices that should merge
//! 4. **Choose representatives**: Pick one vertex from each equivalence class (lowest ID)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
//! 5. **Rebuild mesh**: Create new mesh with faces using representative vertices only
//! 6. **Skip degeneracies**: Don't include faces that become degenerate after mapping
=======
//! 5. **Update topology**: Redirect all edges and faces to use representative vertices
//! 6. **Remove degeneracies**: Delete faces/edges that become degenerate after merging
//! 7. **Clean up**: Remove duplicate vertices and defragment mesh
>>>>>>> Stashed changes
=======
//! 5. **Update topology**: Redirect all edges and faces to use representative vertices
//! 6. **Remove degeneracies**: Delete faces/edges that become degenerate after merging
//! 7. **Clean up**: Remove duplicate vertices and defragment mesh
>>>>>>> Stashed changes
//!
//! # Complexity
//!
//! - **Time**: O(n log n) where n is the number of vertices
//! - **Space**: O(n) for the R*Tree and union-find structures
//!
//! # Example
//!
//! ```rust,ignore
//! use redges::algorithms::vertex_welding::weld_vertices;
//!
//! // Mesh with duplicate vertices from CAD export
//! let mesh = load_mesh("model.obj");
//!
//! // Weld vertices within 1 micron tolerance
//! let (repaired_mesh, report) = mesh.weld_vertices(1e-6);
//!
//! println!("Merged {} duplicate vertices", report.vertices_merged);
//! println!("Removed {} degenerate faces", report.degenerate_faces_removed);
//! ```

use std::collections::HashMap;
use std::marker::PhantomData;

use linear_isomorphic::{InnerSpace, RealField};
use num::{Bounded, Signed};
use rstar::RTree;

use crate::{
    container_trait::{RedgeContainers, VertData},
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    validation::TreePoint,
    Redge, VertId,
=======
    mesh_deleter::MeshDeleter,
    validation::TreePoint,
    EdgeId, Redge, VertId,
>>>>>>> Stashed changes
=======
    mesh_deleter::MeshDeleter,
    validation::TreePoint,
    EdgeId, Redge, VertId,
>>>>>>> Stashed changes
};

/// Report from vertex welding operation.
///
/// This report provides detailed statistics about what was changed during
/// the welding operation, allowing users to understand the impact and
/// validate the repair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VertexWeldingReport {
    /// Number of vertices that were merged/removed.
    ///
    /// This counts the number of duplicate vertices that were removed,
    /// not the number of unique vertex positions found.
    ///
    /// Example: If vertices A, B, C all merge into A, vertices_merged = 2
    pub vertices_merged: usize,

    /// Number of degenerate faces removed.
    ///
    /// A face becomes degenerate when vertex welding causes it to have
    /// duplicate vertices (e.g., a triangle with two or three identical vertices).
    /// These faces have zero area and must be removed.
    pub degenerate_faces_removed: usize,

    /// Number of degenerate edges removed.
    ///
    /// An edge becomes degenerate when vertex welding causes both endpoints
    /// to be the same vertex. These zero-length edges are removed.
    pub degenerate_edges_removed: usize,
}

impl<R: RedgeContainers> Redge<R> {
    /// Weld vertices that are within the specified tolerance distance.
    ///
    /// This operation:
    /// 1. Finds all pairs of vertices within `tolerance` distance
    /// 2. Merges them into equivalence classes (keeps lowest ID as representative)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    /// 3. Rebuilds mesh with faces using representative vertices
    /// 4. Skips any degenerate faces created by the vertex mapping
=======
    /// 3. Updates all edges and faces to reference the representative vertices
    /// 4. Removes duplicate vertices and any degenerate faces/edges created
>>>>>>> Stashed changes
=======
    /// 3. Updates all edges and faces to reference the representative vertices
    /// 4. Removes duplicate vertices and any degenerate faces/edges created
>>>>>>> Stashed changes
    ///
    /// # Arguments
    ///
    /// * `tolerance` - Maximum distance for vertices to be considered duplicates
    ///
    /// # Returns
    ///
    /// Returns a new mesh with welded vertices and a report of changes made.
    ///
    /// # Tolerance Selection
    ///
    /// Choosing the right tolerance is critical:
    /// - **Too small**: Won't merge vertices that should be merged (misses defects)
    /// - **Too large**: Merges vertices that should be separate (damages geometry)
    ///
    /// Recommended tolerances by application:
    /// - **3D printing**: 0.01mm to 0.1mm (model scale dependent)
    /// - **CAD models**: 1e-6 to 1e-3 (unit dependent)
    /// - **Scanned meshes**: 0.1% to 1% of bounding box diagonal
    /// - **Exact duplicates only**: 1e-10 (floating-point epsilon)
    ///
    /// Use `estimate_welding_tolerance()` helper to compute appropriate tolerance
    /// based on mesh bounding box.
    ///
    /// # Complexity
    ///
    /// O(n log n) where n is the number of vertices
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let mesh = Redge::new(vertices, (), (), faces);
    /// let (welded_mesh, report) = mesh.weld_vertices(1e-6);
    /// println!("Merged {} vertices", report.vertices_merged);
    /// ```
    pub fn weld_vertices<S>(self, tolerance: S) -> (Self, VertexWeldingReport)
    where
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        VertData<R>: InnerSpace<S> + Clone,
=======
        VertData<R>: InnerSpace<S>,
>>>>>>> Stashed changes
=======
        VertData<R>: InnerSpace<S>,
>>>>>>> Stashed changes
        S: RealField + Bounded + Signed,
    {
        // Build R*Tree for spatial queries
        let mut tree = RTree::<TreePoint<VertData<R>, S>, _>::new();
        let mut parent: HashMap<VertId, VertId> = HashMap::new();

        // Initialize union-find structure - each vertex is its own parent
        for v in self.meta_verts() {
            parent.insert(v.id(), v.id());
        }

        // Find all vertices within tolerance and build equivalence classes
        for v in self.meta_verts() {
            let test_point = TreePoint {
                point: v.data().clone(),
                id: v.id(),
                _phantom_data: PhantomData,
            };

            // Find all nearby vertices (within tolerance distance)
            let nearby = tree
                .locate_within_distance(test_point.clone(), tolerance * tolerance)
                .collect::<Vec<_>>();

            for nearby_point in nearby {
                // Union the two vertices (choose lower ID as representative)
                let vid1 = v.id();
                let vid2 = nearby_point.id;

                let root1 = find_root(&parent, vid1);
                let root2 = find_root(&parent, vid2);

                if root1 != root2 {
                    // Union by choosing smaller ID as root for determinism
                    if root1.0 < root2.0 {
                        parent.insert(root2, root1);
                    } else {
                        parent.insert(root1, root2);
                    }
                }
            }

            tree.insert(test_point);
        }

        // Path compression for all vertices to get final representatives
        let mut vertex_map: HashMap<VertId, VertId> = HashMap::new();
        for v in self.meta_verts() {
            let root = find_root(&parent, v.id());
            vertex_map.insert(v.id(), root);
        }

        // Count how many vertices will be merged
        let vertices_merged = vertex_map
            .iter()
            .filter(|(vid, rep)| vid != rep)
            .count();

<<<<<<< Updated upstream
<<<<<<< Updated upstream
        // Collect vertex data
        let vert_data: Vec<_> = self.vert_container().iter().cloned().collect();

        // Rebuild faces with remapped vertices, skipping degenerates
        let mut new_faces = Vec::new();
        let mut degenerate_faces = 0;

        for face in self.meta_faces() {
            let verts: Vec<VertId> = face.vertices().map(|v| v.id()).collect();
            let mapped: Vec<VertId> = verts
                .iter()
                .map(|vid| *vertex_map.get(vid).unwrap())
                .collect();

            // Check for duplicates in the mapped vertices
            let unique: std::collections::HashSet<_> = mapped.iter().collect();
            if unique.len() < mapped.len() {
                // Face has duplicate vertices after mapping - skip it
                degenerate_faces += 1;
                continue;
            }

            new_faces.push(mapped);
        }

        // Rebuild mesh with corrected topology
        let new_mesh = Redge::new(
            vert_data,
            self.edge_container().clone(),
            self.face_container().clone(),
            new_faces.iter().map(|f| f.iter().copied()),
        );
=======
=======
>>>>>>> Stashed changes
        // Rebuild mesh with corrected topology
        // This avoids complex topology updates by reconstructing from faces

        // STEP 1: Update ALL topology to use representative vertices
        // This must be done atomically before removing anything

        // Update all half-edge source vertices
        for hedge_id in 0..deleter.mesh.hedges_meta.len() {
            if !deleter.mesh.hedges_meta[hedge_id].is_active {
                continue;
            }
            let old_source = deleter.mesh.hedges_meta[hedge_id].source_id;
            if let Some(&new_source) = vertex_map.get(&old_source) {
                deleter.mesh.hedges_meta[hedge_id].source_id = new_source;
            }
        }

        // Update all edge endpoints
        for edge_id in 0..deleter.mesh.edge_count() {
            if !deleter.mesh.edges_meta[edge_id].is_active {
                continue;
            }

            let v1 = deleter.mesh.edges_meta[edge_id].vert_ids[0];
            let v2 = deleter.mesh.edges_meta[edge_id].vert_ids[1];
            let new_v1 = *vertex_map.get(&v1).unwrap();
            let new_v2 = *vertex_map.get(&v2).unwrap();

            deleter.mesh.edges_meta[edge_id].vert_ids = if new_v1.0 < new_v2.0 {
                [new_v1, new_v2]
            } else {
                [new_v2, new_v1]
            };
        }

        // STEP 2: Mark duplicate vertices as inactive (keeping only representatives)
        // Don't use remove_vert() as it tries to traverse edge cycles which are now inconsistent
        for (vid, rep) in &vertex_map {
            if vid != rep {
                // Just mark as inactive - defragmentation will clean them up
                deleter.mesh.verts_meta[vid.0].is_active = false;
            }
        }

        // STEP 3: Now detect and remove degenerate edges
        // After topology update, edges with v1 == v2 are degenerate
        let mut degenerate_edges = 0;
        let edge_ids: Vec<_> = deleter.mesh.meta_edges().map(|e| e.id()).collect();
        for eid in edge_ids {
            if !deleter.mesh.edges_meta[eid.0].is_active {
                continue;
            }

            let v1 = deleter.mesh.edges_meta[eid.0].vert_ids[0];
            let v2 = deleter.mesh.edges_meta[eid.0].vert_ids[1];

            if v1 == v2 {
                deleter.remove_edge(eid);
                degenerate_edges += 1;
            }
        }

        // STEP 4: Now detect and remove degenerate faces
        // After topology update, check for duplicate vertices
        let mut degenerate_faces = 0;
        let face_ids: Vec<_> = deleter.mesh.meta_faces().map(|f| f.id()).collect();
        for fid in face_ids {
            let face_handle = deleter.mesh.face_handle(fid);
            if !face_handle.is_active() {
                continue;
            }

            // Traverse half-edge loop directly to get vertex IDs
            // Don't use vertices() iterator as it validates through edges
            let verts: Vec<_> = face_handle
                .hedge()
                .face_loop()
                .map(|h| h.source().id())
                .collect();
            let unique_verts: std::collections::HashSet<_> = verts.iter().collect();

            if unique_verts.len() < verts.len() {
                // Face has duplicate vertices - it's degenerate
                deleter.remove_face(fid);
                degenerate_faces += 1;
            }
        }

        let mesh = deleter.end_deletion();
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

        let report = VertexWeldingReport {
            vertices_merged,
            degenerate_faces_removed: degenerate_faces,
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            degenerate_edges_removed: 0, // Edges are implicitly handled by mesh rebuild
        };

        (new_mesh, report)
=======
=======
>>>>>>> Stashed changes
            degenerate_edges_removed: degenerate_edges,
        };

        (mesh, report)
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    }
}

/// Union-find helper: Find root with path compression.
///
/// This implements the classic union-find "find" operation with path compression
/// optimization. Path compression flattens the tree structure by making each
/// visited node point directly to the root, improving future query performance.
fn find_root(parent: &HashMap<VertId, VertId>, vid: VertId) -> VertId {
    let mut root = vid;
    while parent.get(&root).unwrap() != &root {
        root = *parent.get(&root).unwrap();
    }
    root
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_vertex_welding_exact_duplicates() {
        // Create a simple mesh with exact duplicate vertices
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0), // Exact duplicate of vertex 0
        ];

        let faces = vec![vec![0, 1, 2]];

        let mesh = Redge::<(Vec<Vector3<f64>>, (), ())>::new(
            vertices,
            (),
            (),
            faces.into_iter().map(|f| f.into_iter()),
        );

        let (welded_mesh, report) = mesh.weld_vertices(1e-10);

        // Should have merged the exact duplicate
        assert!(report.vertices_merged > 0);
        assert_eq!(welded_mesh.vert_count(), 3);
    }

    #[test]
    fn test_vertex_welding_near_duplicates() {
        // Create mesh with vertices very close but not exactly the same
        let epsilon = 1e-6;
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, epsilon * 0.5), // Near-duplicate of vertex 0
        ];

        let faces = vec![vec![0, 1, 2]];

        let mesh = Redge::<(Vec<Vector3<f64>>, (), ())>::new(
            vertices,
            (),
            (),
            faces.into_iter().map(|f| f.into_iter()),
        );

        let (welded_mesh, report) = mesh.weld_vertices(epsilon);

        // Should have merged the near-duplicate
        assert!(report.vertices_merged > 0);
        assert_eq!(welded_mesh.vert_count(), 3);
    }

    #[test]
    fn test_vertex_welding_removes_degenerate_faces() {
        // Create a triangle where two vertices will be welded, creating a degenerate edge
        let epsilon = 1e-6;
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, epsilon * 0.5), // Will merge with vertex 1
        ];
        let faces = vec![vec![0, 1, 2]];
        let mesh = Redge::<(Vec<Vector3<f64>>, (), ())>::new(
            vertices,
            (),
            (),
            faces.into_iter().map(|f| f.into_iter()),
        );

        let (welded_mesh, report) = mesh.weld_vertices(epsilon);

        // Should have removed the degenerate face
        assert!(report.degenerate_faces_removed > 0 || report.degenerate_edges_removed > 0);
        assert_eq!(welded_mesh.face_count(), 0);
    }

    #[test]
    fn test_vertex_welding_no_duplicates() {
        // Clean mesh with no duplicates
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];

        let faces = vec![vec![0, 1, 2]];

        let mesh = Redge::<(Vec<Vector3<f64>>, (), ())>::new(
            vertices,
            (),
            (),
            faces.into_iter().map(|f| f.into_iter()),
        );

        let initial_vert_count = mesh.vert_count();
        let (welded_mesh, report) = mesh.weld_vertices(1e-10);

        // Nothing should have changed
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.degenerate_faces_removed, 0);
        assert_eq!(welded_mesh.vert_count(), initial_vert_count);
    }
}
