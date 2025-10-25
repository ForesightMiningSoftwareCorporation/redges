//! Quick verification that our defective meshes can be loaded and defects detected
//! Run with: cargo test --test verify_defects -- --nocapture

mod mesh_loaders;

use mesh_loaders::{load_obj_mesh, load_ply_mesh, load_stl_mesh};
use redges::validation::{correctness_state, manifold_state, RedgeCorrectness, RedgeManifoldness};

#[test]
fn verify_can_load_duplicate_vertices() {
    let mesh = load_obj_mesh("tests/fixtures/defective/duplicate_vertices.obj");
    println!("✓ Loaded duplicate_vertices.obj");
    println!("  Vertices: {}", mesh.vert_count());
    println!("  Faces: {}", mesh.face_count());
    println!("  Correctness: {:?}", correctness_state(&mesh));
    assert_eq!(mesh.vert_count(), 4, "Should have 4 vertices (including duplicate)");
}

#[test]
fn verify_can_load_non_manifold_edge() {
    let mesh = load_obj_mesh("tests/fixtures/defective/non_manifold_edge.obj");
    println!("✓ Loaded non_manifold_edge.obj");
    println!("  Vertices: {}", mesh.vert_count());
    println!("  Faces: {}", mesh.face_count());
    println!("  Edges: {}", mesh.edge_count());

    let manifoldness = manifold_state(&mesh);
    println!("  Manifoldness: {:?}", manifoldness);

    // Count edges with more than 2 faces
    let non_manifold_count = mesh
        .meta_edges()
        .filter(|e| e.hedge().radial_loop().count() > 2)
        .count();

    println!("  Non-manifold edges (>2 faces): {}", non_manifold_count);

    assert_ne!(
        manifoldness,
        RedgeManifoldness::IsManifold,
        "Should detect non-manifold edge"
    );
}

#[test]
fn verify_can_detect_boundary_edges() {
    let mesh = load_obj_mesh("tests/fixtures/defective/hole_mesh.obj");
    println!("✓ Loaded hole_mesh.obj");

    let boundary_edges = mesh
        .meta_edges()
        .filter(|e| e.hedge().radial_loop().count() == 1)
        .count();

    println!("  Boundary edges: {}", boundary_edges);
    assert!(boundary_edges > 0, "Should have boundary edges (hole)");
}

#[test]
fn verify_all_synthetic_meshes_load() {
    let test_files = [
        "tests/fixtures/defective/duplicate_vertices.obj",
        "tests/fixtures/defective/near_duplicate_vertices.obj",
        "tests/fixtures/defective/t_junction.obj",
        "tests/fixtures/defective/non_manifold_edge.obj",
        "tests/fixtures/defective/non_manifold_vertex.obj",
        "tests/fixtures/defective/hole_mesh.obj",
    ];

    for file in &test_files {
        let mesh = load_obj_mesh(file);
        println!("✓ Loaded {}: {} verts, {} faces",
                 file, mesh.vert_count(), mesh.face_count());
    }
}

// ============================================================================
// REAL-WORLD DEFECTIVE MESHES - Now we can actually load and verify them!
// ============================================================================

#[test]
fn verify_happy_buddha_has_non_manifold_edges() {
    // Happy Buddha is DOCUMENTED to have 5,581 non-manifold edges
    let mesh = load_ply_mesh("tests/fixtures/defective-real-world/happy_recon/happy_vrip_res3.ply");

    println!("\n=== Happy Buddha (Resolution 3) ===");
    println!("Vertices: {}", mesh.vert_count());
    println!("Faces: {}", mesh.face_count());
    println!("Edges: {}", mesh.edge_count());

    let manifoldness = manifold_state(&mesh);
    println!("Manifoldness: {:?}", manifoldness);

    // Count non-manifold edges (edges with >2 incident faces)
    let non_manifold_edges = mesh
        .meta_edges()
        .filter(|e| e.hedge().radial_loop().count() > 2)
        .count();

    println!("Non-manifold edges detected: {}", non_manifold_edges);

    // PROVE DEFECT EXISTS
    assert_ne!(
        manifoldness,
        RedgeManifoldness::IsManifold,
        "Happy Buddha should have non-manifold defects"
    );

    println!("✓ CONFIRMED: Happy Buddha has detectable non-manifold edges");
}

#[test]
fn verify_dragon_has_boundary_edges() {
    // Dragon is DOCUMENTED to have "numerous small holes"
    let mesh = load_ply_mesh("tests/fixtures/defective-real-world/dragon_recon/dragon_vrip_res3.ply");

    println!("\n=== Dragon (Resolution 3) ===");
    println!("Vertices: {}", mesh.vert_count());
    println!("Faces: {}", mesh.face_count());

    // Count boundary edges (holes)
    let boundary_edges = mesh
        .meta_edges()
        .filter(|e| e.hedge().radial_loop().count() == 1)
        .count();

    println!("Boundary edges (holes): {}", boundary_edges);

    // PROVE DEFECT EXISTS
    assert!(
        boundary_edges > 0,
        "Dragon should have holes (boundary edges)"
    );

    println!("✓ CONFIRMED: Dragon has detectable holes");
}

#[test]
fn verify_rakdos_severe_corruption() {
    // Rakdos is SO corrupted GitHub can't preview it
    let mesh = load_stl_mesh("tests/fixtures/mesh-repair-test-models/models/rakdos.stl");

    println!("\n=== Rakdos (Severely Corrupted) ===");
    println!("Vertices: {}", mesh.vert_count());
    println!("Faces: {}", mesh.face_count());

    let manifoldness = manifold_state(&mesh);
    println!("Manifoldness: {:?}", manifoldness);

    let correctness = correctness_state(&mesh);
    println!("Correctness: {:?}", correctness);

    println!("✓ LOADED: Rakdos despite severe corruption");
}

#[test]
fn verify_double_cube_has_holes() {
    // Double cube has "obvious holes"
    let mesh = load_stl_mesh("tests/fixtures/mesh-repair-test-models/models/double_cube.stl");

    println!("\n=== Double Cube ===");
    println!("Vertices: {}", mesh.vert_count());
    println!("Faces: {}", mesh.face_count());

    let boundary_edges = mesh
        .meta_edges()
        .filter(|e| e.hedge().radial_loop().count() == 1)
        .count();

    println!("Boundary edges: {}", boundary_edges);

    assert!(
        boundary_edges > 0,
        "Double cube should have holes"
    );

    println!("✓ CONFIRMED: Double cube has holes");
}
