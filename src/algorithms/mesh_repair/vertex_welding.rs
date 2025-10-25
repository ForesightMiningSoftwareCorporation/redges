//! Vertex welding algorithm for merging duplicate or nearby vertices.
//!
//! This module implements tolerance-based vertex welding using spatial indexing
//! (R*Tree) and union-find for efficient grouping of duplicate vertices.
//!
//! # Algorithm Overview
//!
//! 1. **Spatial Indexing**: Build R*Tree with all vertex positions for O(log n) queries
//! 2. **Duplicate Detection**: For each vertex, find all neighbors within tolerance
//! 3. **Union-Find Grouping**: Merge vertices into equivalence classes
//! 4. **Representative Selection**: Choose one vertex per group (lowest ID)
//! 5. **Mesh Rebuild**: Reconstruct mesh with faces referencing representatives
//! 6. **Degeneracy Removal**: Skip faces that become degenerate after merging
//!
//! # Complexity
//!
//! - **Time**: O(n log n) where n = number of vertices
//! - **Space**: O(n) for spatial index and union-find structure

use std::collections::HashMap;
use std::marker::PhantomData;

use linear_isomorphic::InnerSpace;
use rstar::RTree;

use crate::{validation::TreePoint, Redge, VertId};

/// Report from vertex welding operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VertexWeldingReport {
    /// Number of vertices that were merged/removed.
    pub vertices_merged: usize,
    /// Number of degenerate faces removed.
    pub degenerate_faces_removed: usize,
}

impl<V> Redge<(Vec<V>, (), ())>
where
    V: Clone + Default + std::fmt::Debug + InnerSpace<f64>,
{
    /// Weld vertices within the specified tolerance distance.
    ///
    /// Merges vertices that are within `tolerance` distance of each other,
    /// using the lowest vertex ID as the representative for each group.
    ///
    /// # Arguments
    ///
    /// * `tolerance` - Maximum distance for vertices to be considered duplicates
    ///
    /// # Returns
    ///
    /// Returns `(new_mesh, report)` where the new mesh has merged vertices
    /// and the report contains statistics about the operation.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let (welded_mesh, report) = mesh.weld_vertices(1e-6);
    /// println!("Merged {} vertices", report.vertices_merged);
    /// ```
    pub fn weld_vertices(self, tolerance: f64) -> (Self, VertexWeldingReport) {
        // Build R*Tree for spatial queries
        let mut tree = RTree::<TreePoint<V, f64>, _>::new();
        let mut parent: HashMap<VertId, VertId> = HashMap::new();

        // Initialize union-find: each vertex starts as its own parent
        for v in self.meta_verts() {
            parent.insert(v.id(), v.id());
        }

        // Find nearby vertices and union them into equivalence classes
        for v in self.meta_verts() {
            let test_point = TreePoint {
                point: v.data().clone(),
                id: v.id(),
                _phantom_data: PhantomData,
            };

            // Query R*Tree for all vertices within tolerance distance
            let nearby = tree
                .locate_within_distance(test_point.clone(), tolerance * tolerance)
                .collect::<Vec<_>>();

            for nearby_point in nearby {
                let vid1 = v.id();
                let vid2 = nearby_point.id;

                let root1 = find_root(&parent, vid1);
                let root2 = find_root(&parent, vid2);

                // Union by choosing smaller ID as root (for determinism)
                if root1 != root2 {
                    if root1.0 < root2.0 {
                        parent.insert(root2, root1);
                    } else {
                        parent.insert(root1, root2);
                    }
                }
            }

            tree.insert(test_point);
        }

        // Build final vertex mapping with path compression
        let mut vertex_map: HashMap<VertId, VertId> = HashMap::new();
        for v in self.meta_verts() {
            let root = find_root(&parent, v.id());
            vertex_map.insert(v.id(), root);
        }

        let vertices_merged = vertex_map.iter().filter(|(vid, rep)| vid != rep).count();

        // Extract mesh data
        let (old_verts, face_indices, _face_data) = self.to_face_list();

        // Build mapping from old vertex ID to new compact vertex ID
        // Only include representative vertices
        let mut old_to_new: HashMap<VertId, usize> = HashMap::new();
        let mut new_verts = Vec::new();

        for vid in vertex_map.values() {
            if !old_to_new.contains_key(vid) {
                old_to_new.insert(*vid, new_verts.len());
                new_verts.push(old_verts[vid.0].clone());
            }
        }

        // Rebuild faces with remapped vertices, skipping degenerates
        let mut new_faces = Vec::new();
        let mut degenerate_faces = 0;

        for face_verts in face_indices {
            // Map face vertices to their representatives, then to new compact IDs
            let mapped: Vec<usize> = face_verts
                .iter()
                .map(|&idx| {
                    let rep = vertex_map.get(&VertId(idx)).unwrap();
                    *old_to_new.get(rep).unwrap()
                })
                .collect();

            // Check for duplicate vertices (degenerate face)
            let unique: std::collections::HashSet<_> = mapped.iter().collect();
            if unique.len() < mapped.len() {
                degenerate_faces += 1;
                continue;
            }

            new_faces.push(mapped);
        }

        // Rebuild mesh with welded vertices
        let new_mesh = Redge::new(
            new_verts,
            (),
            (),
            new_faces.iter().map(|f| f.iter().copied()),
        );

        let report = VertexWeldingReport {
            vertices_merged,
            degenerate_faces_removed: degenerate_faces,
        };

        (new_mesh, report)
    }
}

/// Find root of vertex in union-find structure with path compression.
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

    type TestMesh = Redge<(Vec<Vector3<f64>>, (), ())>;

    #[test]
    fn test_weld_exact_duplicates() {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0), // Exact duplicate of vertex 0
        ];
        let faces = vec![vec![0, 1, 2]];
        let mesh: TestMesh = Redge::new(vertices, (), (), faces.into_iter().map(|f| f.into_iter()));

        let (welded, report) = mesh.weld_vertices(1e-10);

        assert_eq!(report.vertices_merged, 1);
        assert_eq!(welded.vert_count(), 3);
    }

    #[test]
    fn test_weld_near_duplicates() {
        let epsilon = 1e-6;
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, epsilon * 0.5), // Near-duplicate
        ];
        let faces = vec![vec![0, 1, 2]];
        let mesh: TestMesh = Redge::new(vertices, (), (), faces.into_iter().map(|f| f.into_iter()));

        let (welded, report) = mesh.weld_vertices(epsilon);

        assert_eq!(report.vertices_merged, 1);
        assert_eq!(welded.vert_count(), 3);
    }

    #[test]
    fn test_weld_removes_degenerate_faces() {
        let epsilon = 1e-6;
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, epsilon * 0.5), // Will merge with vertex 1
        ];
        let faces = vec![vec![0, 1, 2]];
        let mesh: TestMesh = Redge::new(vertices, (), (), faces.into_iter().map(|f| f.into_iter()));

        let (welded, report) = mesh.weld_vertices(epsilon);

        assert!(report.degenerate_faces_removed > 0);
        assert_eq!(welded.face_count(), 0);
    }

    #[test]
    fn test_weld_no_duplicates() {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let faces = vec![vec![0, 1, 2]];
        let mesh: TestMesh = Redge::new(vertices, (), (), faces.into_iter().map(|f| f.into_iter()));

        let (welded, report) = mesh.weld_vertices(1e-10);

        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.degenerate_faces_removed, 0);
        assert_eq!(welded.vert_count(), 3);
    }
}
