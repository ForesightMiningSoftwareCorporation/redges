//! Alternative vertex welding implementation using mesh rebuild approach

use std::collections::HashMap;
use std::marker::PhantomData;

use linear_isomorphic::{InnerSpace, RealField};
use num::{Bounded, Signed};
use rstar::RTree;

use crate::{
    container_trait::{RedgeContainers, VertData},
    validation::TreePoint,
    Redge, VertId,
};

/// Simplified report from vertex welding
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeldReport {
    pub vertices_merged: usize,
    pub degenerate_faces_removed: usize,
}

impl<R: RedgeContainers> Redge<R> {
    /// Weld vertices by rebuilding the mesh with updated face indices
    pub fn weld_vertices_v2<S>(self, tolerance: S) -> (Self, WeldReport)
    where
        VertData<R>: InnerSpace<S> + Clone,
        S: RealField + Bounded + Signed,
    {
        // Build R*Tree and find duplicate vertices
        let mut tree = RTree::<TreePoint<VertData<R>, S>, _>::new();
        let mut parent: HashMap<VertId, VertId> = HashMap::new();

        // Initialize union-find
        for v in self.meta_verts() {
            parent.insert(v.id(), v.id());
        }

        // Find nearby vertices and union them
        for v in self.meta_verts() {
            let test_point = TreePoint {
                point: v.data().clone(),
                id: v.id(),
                _phantom_data: PhantomData,
            };

            let nearby = tree
                .locate_within_distance(test_point.clone(), tolerance * tolerance)
                .collect::<Vec<_>>();

            for nearby_point in nearby {
                let vid1 = v.id();
                let vid2 = nearby_point.id;

                let root1 = find_root(&parent, vid1);
                let root2 = find_root(&parent, vid2);

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

        // Get final vertex mapping
        let mut vertex_map: HashMap<VertId, VertId> = HashMap::new();
        for v in self.meta_verts() {
            let root = find_root(&parent, v.id());
            vertex_map.insert(v.id(), root);
        }

        let vertices_merged = vertex_map.iter().filter(|(vid, rep)| vid != rep).count();

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

            // Check for duplicates
            let unique: std::collections::HashSet<_> = mapped.iter().collect();
            if unique.len() < mapped.len() {
                degenerate_faces += 1;
                continue;
            }

            new_faces.push(mapped);
        }

        // Rebuild mesh
        let new_mesh = Redge::new(
            vert_data,
            self.edge_container().clone(),
            self.face_container().clone(),
            new_faces.iter().map(|f| f.iter().copied()),
        );

        let report = WeldReport {
            vertices_merged,
            degenerate_faces_removed: degenerate_faces,
        };

        (new_mesh, report)
    }
}

fn find_root(parent: &HashMap<VertId, VertId>, mut vid: VertId) -> VertId {
    let mut root = vid;
    while parent.get(&root).unwrap() != &root {
        root = *parent.get(&root).unwrap();
    }
    root
}
