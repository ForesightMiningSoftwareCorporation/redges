//! Tests that demonstrate real mesh defects exist and require repair.
//!
//! **Purpose**: Prove that mesh repair is NECESSARY for real-world geometry.
//!
//! This test suite:
//! 1. Loads meshes with **documented, known defects** from academic sources
//! 2. **Proves these defects exist** by detecting them programmatically
//! 3. Shows that repair operations **fix these defects**
//!
//! We use meshes that are explicitly broken and well-cited in research:
//! - Happy Buddha: 153 self-intersections, 5581 non-manifold edges (Stanford)
//! - Dragon: Numerous holes from scanning (Stanford)
//! - Rakdos: Severely corrupted with many non-manifold holes
//! - Synthetic defects: Controlled test cases for specific issues

use nalgebra::Vector3;
use redges::{
    algorithms::mesh_repair::*,
    validation::{correctness_state, manifold_state, RedgeCorrectness, RedgeManifoldness},
    wavefront_loader::ObjData,
    Redge,
};

type TestMesh = Redge<(Vec<Vector3<f64>>, (), ())>;

/// Helper function to load OBJ file into TestMesh
fn load_obj_mesh(path: &str) -> TestMesh {
    let obj_data = ObjData::from_disk_file(path);
    let vertices: Vec<_> = obj_data.vertices.iter().map(|v| v.map(|s| s as f64)).collect();
    let indices: Vec<_> = obj_data
        .vertex_face_indices
        .iter()
        .map(|l| l.iter().map(|&i| i as usize).collect::<Vec<_>>())
        .collect();

    Redge::new(
        vertices,
        (),
        (),
        indices.iter().map(|f| f.iter().copied()),
    )
}

// ============================================================================
// SYNTHETIC DEFECTS - Controlled test cases
// ============================================================================

#[test]
fn test_exact_duplicate_vertices_defect_exists() {
    let obj_path = "tests/fixtures/defective/duplicate_vertices.obj";
    let mesh: TestMesh = load_obj_mesh(obj_path);

    // PROVE DEFECT EXISTS: Should have 4 vertices including the duplicate
    assert_eq!(
        mesh.vert_count(),
        4,
        "DEFECT: Mesh has duplicate vertices (4 instead of 3 unique)"
    );

    // NOW FIX IT
    let (welded_mesh, report) = mesh.weld_vertices(1e-10);

    // PROVE FIX WORKED
    assert!(
        report.vertices_merged > 0,
        "Repair should have merged duplicate vertices"
    );
    assert_eq!(
        welded_mesh.vert_count(),
        3,
        "After repair: 3 unique vertices"
    );

    // Verify mesh is still correct after repair
    assert_eq!(
        correctness_state(&welded_mesh),
        RedgeCorrectness::Correct,
        "Mesh should be correct after repair"
    );

    println!("✓ Exact duplicate defect detected and fixed");
    println!("  Before: {} vertices", 4);
    println!("  After:  {} vertices", welded_mesh.vert_count());
    println!("  Merged: {}", report.vertices_merged);
}

#[test]
fn test_near_duplicate_vertices_within_tolerance() {
    let obj_path = "tests/fixtures/defective/near_duplicate_vertices.obj";
    let mesh: TestMesh = load_obj_mesh(obj_path);

    // PROVE DEFECT EXISTS
    assert_eq!(mesh.vert_count(), 4, "DEFECT: Has near-duplicate vertices");

    // Very tight tolerance should NOT fix it (proves vertices aren't exactly the same)
    let (tight_welded, tight_report) = mesh.clone().weld_vertices(1e-10);
    assert_eq!(
        tight_report.vertices_merged, 0,
        "Vertices are NOT exact duplicates"
    );
    assert_eq!(tight_welded.vert_count(), 4);

    // Reasonable tolerance SHOULD fix it (proves they're within tolerance)
    let (welded_mesh, report) = mesh.weld_vertices(0.001);
    assert!(
        report.vertices_merged > 0,
        "REPAIR: Vertices within 0.001 should be merged"
    );
    assert_eq!(welded_mesh.vert_count(), 3);

    println!("✓ Near-duplicate defect (tolerance-based) detected and fixed");
    println!("  Tight tolerance (1e-10): {} merged", tight_report.vertices_merged);
    println!("  Loose tolerance (1e-3):  {} merged", report.vertices_merged);
}

#[test]
fn test_t_junction_defect_common_in_cad() {
    // T-junctions are a VERY common defect in CAD exports
    let obj_path = "tests/fixtures/defective/t_junction.obj";
    let mesh: TestMesh = load_obj_mesh(obj_path);

    let original_vert_count = mesh.vert_count();

    // PROVE DEFECT EXISTS: There should be near-duplicate vertices at the junction
    let (welded_mesh, report) = mesh.weld_vertices(0.0001);

    assert!(
        report.vertices_merged > 0,
        "DEFECT: T-junction has near-duplicate vertices"
    );
    assert!(
        welded_mesh.vert_count() < original_vert_count,
        "REPAIR: Welding should reduce vertex count"
    );

    println!("✓ T-junction defect detected and fixed");
    println!("  Before: {} vertices", original_vert_count);
    println!("  After:  {} vertices", welded_mesh.vert_count());
    println!("  This is a REAL problem in CAD/boolean operations");
}

#[test]
fn test_non_manifold_edge_violates_2_manifold_property() {
    // This mesh has an edge shared by 3 faces - violates 2-manifold
    let obj_path = "tests/fixtures/defective/non_manifold_edge.obj";
    let mesh: TestMesh = load_obj_mesh(obj_path);

    // PROVE DEFECT EXISTS
    let manifoldness_before = manifold_state(&mesh);

    // Check if it's non-manifold (should not be IsManifold)
    assert_ne!(
        manifoldness_before,
        RedgeManifoldness::IsManifold,
        "DEFECT: Mesh should not be manifold (has edge with 3+ faces)"
    );

    // Count the actual non-manifold edges
    let non_manifold_edges = mesh
        .meta_edges()
        .filter(|e| e.hedge().radial_loop().count() > 2)
        .count();

    assert!(
        non_manifold_edges > 0,
        "Should detect at least one non-manifold edge"
    );

    println!("✓ Non-manifold edge defect detected");
    println!("  Non-manifold edges: {}", non_manifold_edges);
    println!("  Manifoldness: {:?}", manifoldness_before);
    println!("  TODO: Implement repair for this defect");

    // TODO: When non-manifold repair is implemented:
    // let mut deleter = MeshDeleter::start_deletion(mesh);
    // let report = deleter.repair_nonmanifold_topology();
    // let repaired_mesh = deleter.end_deletion();
    // assert_eq!(repaired_mesh.is_manifold(), RedgeManifoldness::IsManifold);
}

#[test]
fn test_non_manifold_vertex_bowtie_configuration() {
    // Bowtie: vertex shared by two disconnected triangle fans
    let obj_path = "tests/fixtures/defective/non_manifold_vertex.obj";
    let mesh: TestMesh = load_obj_mesh(obj_path);

    // PROVE DEFECT EXISTS
    let manifoldness = manifold_state(&mesh);

    assert_ne!(
        manifoldness,
        RedgeManifoldness::IsManifold,
        "DEFECT: Bowtie configuration should NOT be manifold"
    );

    println!("✓ Non-manifold vertex (bowtie) defect detected");
    println!("  Manifoldness: {:?}", manifoldness);
    println!("  This is a topological defect requiring vertex splitting");
}

#[test]
fn test_hole_creates_boundary_edges() {
    // Mesh with missing faces = boundary edges = hole
    let obj_path = "tests/fixtures/defective/hole_mesh.obj";
    let mesh: TestMesh = load_obj_mesh(obj_path);

    // PROVE DEFECT EXISTS: Count boundary edges (edges with only 1 face)
    let boundary_edge_count = mesh
        .meta_edges()
        .filter(|e| {
            let radial_count = e.hedge().radial_loop().count();
            radial_count == 1
        })
        .count();

    assert!(
        boundary_edge_count > 0,
        "DEFECT: Mesh has hole (boundary edges exist)"
    );

    println!("✓ Hole defect detected");
    println!("  Boundary edges: {}", boundary_edge_count);
    println!("  Mesh is NOT watertight/solid");
    println!("  TODO: Implement hole filling");
}

// ============================================================================
// REAL-WORLD DEFECTS - Documented broken meshes from research
// ============================================================================

// NOTE: These tests are marked #[ignore] because they require PLY format support
// or conversion to OBJ. They document what SHOULD be tested once we have PLY loading.

#[test]
#[ignore = "Requires PLY loader - Happy Buddha is in PLY format"]
fn test_happy_buddha_has_documented_defects() {
    // Happy Buddha is FAMOUS for having:
    // - 153 self-intersections
    // - 5,581 non-manifold edges
    // This is cited in academic mesh repair papers!

    // TODO: Load PLY file
    // let ply_path = "tests/fixtures/defective-real-world/happy_recon/happy_vrip_res2.ply";
    // let mesh = load_ply_mesh(ply_path);

    println!("Happy Buddha (Stanford) - KNOWN DEFECTS:");
    println!("  • 153 self-intersections (documented)");
    println!("  • 5,581 non-manifold edges (documented)");
    println!("  • Small bridges from space carving");
    println!("  • Topological genus larger than appears");
    println!("  Source: Stanford 3D Scanning Repository");
    println!("  Citation: Used in mesh repair benchmarks");

    // Once PLY loading works:
    // let non_manifold_edges = mesh.meta_edges()
    //     .filter(|e| e.hedge().radial_loop().count() > 2)
    //     .count();
    // assert!(non_manifold_edges > 5000, "Should have ~5581 non-manifold edges");
}

#[test]
#[ignore = "Requires PLY loader - Dragon is in PLY format"]
fn test_dragon_has_numerous_holes() {
    // Stanford Dragon is documented to have "numerous small holes"
    // from scanning artifacts

    // TODO: Load PLY file
    // let ply_path = "tests/fixtures/defective-real-world/dragon_recon/dragon_vrip_res2.ply";
    // let mesh = load_ply_mesh(ply_path);

    println!("Dragon (Stanford) - KNOWN DEFECTS:");
    println!("  • Numerous small holes (documented)");
    println!("  • Incomplete surface from scanning");
    println!("  • 566,098 vertices, 1,132,830 triangles");
    println!("  Source: Stanford 3D Scanning Repository");

    // Once PLY loading works:
    // let boundary_edges = mesh.meta_edges()
    //     .filter(|e| e.hedge().radial_loop().count() == 1)
    //     .count();
    // assert!(boundary_edges > 0, "Should have boundary edges (holes)");
}

#[test]
#[ignore = "Requires STL loader"]
fn test_rakdos_severely_corrupted() {
    // Rakdos is SO corrupted that GitHub 3D preview can't even render it
    // Documented as having "many tiny non-manifold holes"

    // TODO: Load STL file
    // let stl_path = "tests/fixtures/mesh-repair-test-models/models/rakdos.stl";
    // let mesh = load_stl_mesh(stl_path);

    println!("Rakdos - KNOWN DEFECTS:");
    println!("  • Many tiny non-manifold holes (documented)");
    println!("  • SO corrupted GitHub can't preview it");
    println!("  • ~10MB file size");
    println!("  Source: caretdashcaret/MeshRepairTestModels");
    println!("  This is a STRESS TEST for extreme corruption");
}

#[test]
#[ignore = "Requires STL loader"]
fn test_double_cube_has_obvious_holes() {
    // Simple case: documented to have "obvious holes that need to be filled"

    // TODO: Load STL file
    // let stl_path = "tests/fixtures/mesh-repair-test-models/models/double_cube.stl";
    // let mesh = load_stl_mesh(stl_path);

    println!("Double Cube - KNOWN DEFECTS:");
    println!("  • Obvious holes that need filling (documented)");
    println!("  • Simple test case for hole detection");
    println!("  Source: caretdashcaret/MeshRepairTestModels");
}

// ============================================================================
// REPAIR VALIDATION - Prove fixes work
// ============================================================================

#[test]
fn test_vertex_welding_report_accuracy() {
    // Verify that repair reports are accurate
    let obj_path = "tests/fixtures/defective/duplicate_vertices.obj";
    let mesh: TestMesh = load_obj_mesh(obj_path);

    let original_vert_count = mesh.vert_count();
    let (welded_mesh, report) = mesh.weld_vertices(1e-10);

    // Report should match actual changes
    let actual_removed = original_vert_count - welded_mesh.vert_count();
    assert_eq!(
        report.vertices_merged as usize,
        actual_removed,
        "Report should accurately reflect vertices merged"
    );

    // Mesh should still be valid after repair
    let correctness = correctness_state(&welded_mesh);
    assert_eq!(
        correctness,
        RedgeCorrectness::Correct,
        "Repaired mesh should be topologically correct"
    );

    println!("✓ Repair report validation passed");
    println!("  Reported merged: {}", report.vertices_merged);
    println!("  Actual removed:  {}", actual_removed);
    println!("  Final mesh correctness: {:?}", correctness);
}

#[test]
fn test_tolerance_sensitivity() {
    // Show that tolerance matters for near-duplicates
    let obj_path = "tests/fixtures/defective/near_duplicate_vertices.obj";
    let mesh: TestMesh = load_obj_mesh(obj_path);

    println!("Testing tolerance sensitivity:");

    for tolerance in [1e-10, 1e-6, 1e-4, 1e-3, 1e-2] {
        let (_, report) = mesh.clone().weld_vertices(tolerance);
        println!("  Tolerance {:.0e}: {} vertices merged",
                 tolerance, report.vertices_merged);
    }

    // Very tight: nothing merged
    let (_, tight_report) = mesh.clone().weld_vertices(1e-10);
    assert_eq!(tight_report.vertices_merged, 0);

    // Reasonable: should merge near-duplicates
    let (_, loose_report) = mesh.weld_vertices(1e-3);
    assert!(loose_report.vertices_merged > 0);

    println!("✓ Tolerance affects vertex welding as expected");
}

// ============================================================================
// DOCUMENTATION TESTS - Explain WHY these defects matter
// ============================================================================

#[test]
fn document_why_duplicate_vertices_are_bad() {
    println!("\n=== Why Duplicate Vertices Are a Problem ===");
    println!("• Wastes memory (duplicate data)");
    println!("• Breaks watertightness (creates T-junctions)");
    println!("• Causes gaps in 3D printing");
    println!("• Fails boolean operations");
    println!("• Prevents proper normal calculation");
    println!("• Common in: CAD exports, merging meshes, precision errors\n");
}

#[test]
fn document_why_non_manifold_is_bad() {
    println!("\n=== Why Non-Manifold Geometry Is a Problem ===");
    println!("• Violates 2-manifold property (required for many algorithms)");
    println!("• Cannot determine inside/outside (no solid volume)");
    println!("• 3D printing fails (slicer can't process)");
    println!("• Simulation fails (FEM requires manifold)");
    println!("• Subdivision surfaces break");
    println!("• Common in: Boolean ops, mesh merging, scanning errors\n");
}

#[test]
fn document_why_holes_are_bad() {
    println!("\n=== Why Holes Are a Problem ===");
    println!("• Mesh is not watertight/solid");
    println!("• Cannot 3D print (no volume)");
    println!("• Cannot compute volume/mass");
    println!("• Simulation fails");
    println!("• Texture mapping breaks at boundaries");
    println!("• Common in: 3D scanning, incomplete reconstruction\n");
}
