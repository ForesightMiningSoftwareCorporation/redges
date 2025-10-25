//! Comprehensive test suite for mesh repair operations.
//!
//! This test suite covers:
//! - Vertex welding on synthetic and real-world meshes
//! - Non-manifold edge and vertex repair
//! - Integration with real-world test datasets

use nalgebra::Vector3;
use redges::Redge;

type TestMesh = Redge<(Vec<Vector3<f64>>, (), ())>;

/// Helper to create a simple triangle mesh
fn create_triangle() -> TestMesh {
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.5, 1.0, 0.0),
    ];
    let faces = vec![vec![0, 1, 2]];
    Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    )
}

/// Helper to create a cube with 8 vertices and 12 triangular faces
fn create_cube() -> TestMesh {
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0), // 0
        Vector3::new(1.0, 0.0, 0.0), // 1
        Vector3::new(1.0, 1.0, 0.0), // 2
        Vector3::new(0.0, 1.0, 0.0), // 3
        Vector3::new(0.0, 0.0, 1.0), // 4
        Vector3::new(1.0, 0.0, 1.0), // 5
        Vector3::new(1.0, 1.0, 1.0), // 6
        Vector3::new(0.0, 1.0, 1.0), // 7
    ];
    let faces = vec![
        // Bottom face (z=0)
        vec![0, 2, 1],
        vec![0, 3, 2],
        // Top face (z=1)
        vec![4, 5, 6],
        vec![4, 6, 7],
        // Front face (y=0)
        vec![0, 1, 5],
        vec![0, 5, 4],
        // Back face (y=1)
        vec![2, 3, 7],
        vec![2, 7, 6],
        // Left face (x=0)
        vec![0, 4, 7],
        vec![0, 7, 3],
        // Right face (x=1)
        vec![1, 2, 6],
        vec![1, 6, 5],
    ];
    Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    )
}

// ============================================================================
// Vertex Welding Tests
// ============================================================================

#[test]
fn test_vertex_welding_exact_duplicates() {
    // Create mesh with exact duplicate vertices
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.5, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 0.0), // Exact duplicate of vertex 0
    ];
    let faces = vec![vec![0, 1, 2], vec![1, 3, 2]];
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, report) = mesh.weld_vertices(1e-10);

    assert!(report.vertices_merged > 0);
    assert_eq!(welded.vert_count(), 3); // Should have only 3 unique vertices
}

#[test]
fn test_vertex_welding_near_duplicates() {
    // Create mesh with vertices very close but not exactly the same
    let epsilon = 1e-6;
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.5, 1.0, 0.0),
        Vector3::new(0.0, 0.0, epsilon * 0.5), // Near-duplicate of vertex 0
    ];
    let faces = vec![vec![0, 1, 2]];
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, report) = mesh.weld_vertices(epsilon);

    assert!(report.vertices_merged > 0);
    assert_eq!(welded.vert_count(), 3);
}

#[test]
fn test_vertex_welding_no_duplicates() {
    let mesh = create_cube();
    let initial_vert_count = mesh.vert_count();

    let (welded, report) = mesh.weld_vertices(1e-10);

    assert_eq!(report.vertices_merged, 0);
    assert_eq!(welded.vert_count(), initial_vert_count);
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
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, report) = mesh.weld_vertices(epsilon);

    assert!(report.degenerate_faces_removed > 0 || report.degenerate_edges_removed > 0);
}

#[test]
fn test_vertex_welding_tolerance_sensitivity() {
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.01, 0.0), // 0.01 units away from vertex 0
    ];
    let faces = vec![vec![0, 1, 2]];

    // Test with small tolerance - should not merge
    let mesh1 = Redge::new(
        vertices.clone(),
        (),
        (),
        faces.clone().into_iter().map(|f| f.into_iter()),
    );
    let (_, report1) = mesh1.weld_vertices(0.001);
    assert_eq!(report1.vertices_merged, 0);

    // Test with large tolerance - should merge
    let mesh2 = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );
    let (_, report2) = mesh2.weld_vertices(0.02);
    assert!(report2.vertices_merged > 0);
}

#[test]
fn test_vertex_welding_preserves_clean_geometry() {
    let mesh = create_cube();
    let original_face_count = mesh.face_count();

    let (welded, report) = mesh.weld_vertices(1e-10);

    assert_eq!(report.vertices_merged, 0);
    assert_eq!(report.degenerate_faces_removed, 0);
    assert_eq!(welded.face_count(), original_face_count);
}

#[test]
fn test_vertex_welding_multiple_duplicates() {
    // Create mesh with multiple sets of duplicate vertices
    let epsilon = 1e-6;
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, epsilon * 0.5), // Dup of 0
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, epsilon * 0.5), // Dup of 2
        Vector3::new(0.5, 1.0, 0.0),
        Vector3::new(0.5, 1.0, epsilon * 0.5), // Dup of 4
    ];
    let faces = vec![vec![0, 2, 4], vec![1, 3, 5]];
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, report) = mesh.weld_vertices(epsilon);

    assert!(report.vertices_merged >= 3); // At least 3 merges
    assert_eq!(welded.vert_count(), 3); // Should end up with 3 unique vertices
}

#[test]
fn test_vertex_welding_chain_merging() {
    // Test case where A ≈ B and B ≈ C (transitivity)
    let epsilon = 1e-6;
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(epsilon * 0.8, 0.0, 0.0), // Close to 0
        Vector3::new(epsilon * 1.6, 0.0, 0.0), // Close to 1, farther from 0
        Vector3::new(1.0, 0.0, 0.0),
    ];
    let faces = vec![vec![0, 1, 3], vec![1, 2, 3]];
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, _report) = mesh.weld_vertices(epsilon);

    // Vertices 0, 1, 2 should all merge into one
    assert_eq!(welded.vert_count(), 2);
}

#[test]
fn test_vertex_welding_boundary_vertices() {
    // Create a mesh with a hole (boundary vertices)
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    ];
    // Only two triangles, creating a boundary
    let faces = vec![vec![0, 1, 2], vec![0, 2, 3]];
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, report) = mesh.weld_vertices(1e-10);

    assert_eq!(report.vertices_merged, 0);
    assert_eq!(welded.vert_count(), 4);
}

// ============================================================================
// Non-Manifold Repair Tests
// ============================================================================

#[test]
fn test_split_manifold_edge_noop() {
    use redges::mesh_deleter::MeshDeleter;

    let mesh = create_cube();
    let mut deleter = MeshDeleter::start_deletion(mesh);

    // Pick any edge - in a cube all edges are manifold (2 incident faces)
    let edge_id = deleter.mesh.meta_edges().next().unwrap().id();

    let new_edges = deleter.split_nonmanifold_edge(edge_id);

    assert_eq!(new_edges.len(), 0); // No new edges created for manifold edge
}

#[test]
fn test_repair_manifold_mesh_noop() {
    use redges::mesh_deleter::MeshDeleter;

    let mesh = create_cube();
    let mut deleter = MeshDeleter::start_deletion(mesh);

    let report = deleter.repair_nonmanifold_topology();

    assert_eq!(report.edges_split, 0);
    assert_eq!(report.vertices_split, 0);
}

// ============================================================================
// Integration Tests with Sample Datasets
// ============================================================================

#[test]
#[ignore] // Only run when test data is available
fn test_stanford_bunny_vertex_welding() {
    use redges::wavefront_loader::ObjData;
    use std::path::Path;

    let path = Path::new("tests/fixtures/stanford/bunny.obj");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }

    let obj_data = ObjData::from_disk_file(path.to_str().unwrap());
    let vertices: Vec<_> = obj_data.vertices.iter().map(|v| v.map(|s| s as f64)).collect();
    let indices: Vec<_> = obj_data
        .vertex_face_indices
        .iter()
        .map(|l| l.iter().map(|&i| i as usize).collect::<Vec<_>>())
        .collect();

    let mesh = Redge::<(_, (), ())>::new(
        vertices,
        (),
        (),
        indices.iter().map(|f| f.iter().copied()),
    );

    let initial_vert_count = mesh.vert_count();
    let (welded, report) = mesh.weld_vertices(1e-6);

    println!(
        "Stanford Bunny: {} vertices -> {} vertices ({} merged)",
        initial_vert_count,
        welded.vert_count(),
        report.vertices_merged
    );

    assert!(welded.vert_count() <= initial_vert_count);
}

#[test]
#[ignore]
fn test_cow_vertex_welding() {
    use redges::wavefront_loader::ObjData;
    use std::path::Path;

    let path = Path::new("tests/fixtures/stanford/cow.obj");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }

    let obj_data = ObjData::from_disk_file(path.to_str().unwrap());
    let vertices: Vec<_> = obj_data.vertices.iter().map(|v| v.map(|s| s as f64)).collect();
    let indices: Vec<_> = obj_data
        .vertex_face_indices
        .iter()
        .map(|l| l.iter().map(|&i| i as usize).collect::<Vec<_>>())
        .collect();

    let mesh = Redge::<(_, (), ())>::new(
        vertices,
        (),
        (),
        indices.iter().map(|f| f.iter().copied()),
    );

    let (_, report) = mesh.weld_vertices(1e-6);

    println!("Cow model repair report: {:?}", report);
}

// ============================================================================
// Stress Tests
// ============================================================================

#[test]
fn test_vertex_welding_large_mesh() {
    // Create a larger mesh (100x100 grid = 10000 vertices)
    let grid_size = 100;
    let mut vertices = Vec::new();
    let mut faces = Vec::new();

    for i in 0..grid_size {
        for j in 0..grid_size {
            vertices.push(Vector3::new(i as f64, j as f64, 0.0));
        }
    }

    for i in 0..grid_size - 1 {
        for j in 0..grid_size - 1 {
            let idx = i * grid_size + j;
            faces.push(vec![idx, idx + 1, idx + grid_size + 1]);
            faces.push(vec![idx, idx + grid_size + 1, idx + grid_size]);
        }
    }

    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, report) = mesh.weld_vertices(1e-10);

    assert_eq!(report.vertices_merged, 0); // No duplicates in grid
    assert_eq!(welded.vert_count(), grid_size * grid_size);
}

#[test]
fn test_vertex_welding_all_duplicates() {
    // Extreme case: all vertices are at the same location
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 0.0),
    ];
    let faces = vec![vec![0, 1, 2], vec![1, 2, 3]];
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, report) = mesh.weld_vertices(1e-10);

    assert_eq!(welded.vert_count(), 1); // All merge into one
    assert_eq!(report.vertices_merged, 3);
    // All faces should be degenerate
    assert_eq!(welded.face_count(), 0);
}

#[test]
fn test_vertex_welding_empty_mesh() {
    let vertices: Vec<Vector3<f64>> = vec![];
    let faces: Vec<Vec<usize>> = vec![];
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, report) = mesh.weld_vertices(1e-10);

    assert_eq!(report.vertices_merged, 0);
    assert_eq!(welded.vert_count(), 0);
    assert_eq!(welded.face_count(), 0);
}

#[test]
fn test_vertex_welding_single_triangle() {
    let mesh = create_triangle();
    let (welded, report) = mesh.weld_vertices(1e-10);

    assert_eq!(report.vertices_merged, 0);
    assert_eq!(welded.vert_count(), 3);
    assert_eq!(welded.face_count(), 1);
}

// ============================================================================
// Correctness Validation Tests
// ============================================================================

#[test]
fn test_welded_mesh_is_valid() {
    use redges::validation::{correctness_state, RedgeCorrectness};

    let epsilon = 1e-6;
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.5, 1.0, 0.0),
        Vector3::new(0.0, 0.0, epsilon * 0.5), // Duplicate
    ];
    let faces = vec![vec![0, 1, 2]];
    let mesh = Redge::new(
        vertices,
        (),
        (),
        faces.into_iter().map(|f| f.into_iter()),
    );

    let (welded, _) = mesh.weld_vertices(epsilon);

    // Verify the welded mesh passes correctness checks
    assert_eq!(correctness_state(&welded), RedgeCorrectness::Correct);
}

#[test]
fn test_vertex_welding_preserves_manifoldness() {
    use redges::validation::{manifold_state, RedgeManifoldness};

    let mesh = create_cube();

    // Verify original is manifold
    assert_eq!(manifold_state(&mesh), RedgeManifoldness::IsManifold);

    let (welded, _) = mesh.weld_vertices(1e-10);

    // Verify welded mesh is still manifold
    assert_eq!(manifold_state(&welded), RedgeManifoldness::IsManifold);
}
