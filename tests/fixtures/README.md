# Test Mesh Fixtures - Defective Meshes

This directory contains **defective meshes with known problems** used for validating mesh repair algorithms.

The goal is to demonstrate that mesh repair is NECESSARY - these are not clean, perfect meshes. They are real-world broken geometry that needs fixing.

## Directory Structure

```
fixtures/
├── defective/                  # Synthetic meshes with specific defects
├── defective-real-world/       # Real scanned meshes with documented defects
└── mesh-repair-test-models/    # Known broken meshes from mesh repair research
```

## Real-World Defective Meshes (`defective-real-world/`)

These are **actual scanned meshes from Stanford 3D Scanning Repository that are explicitly documented as having defects**.

### Happy Buddha (`happy_recon/`)
- **Source**: Stanford 3D Scanning Repository
- **Known Defects**: 
  - Contains **153 self-intersections**
  - Contains **5,581 non-manifold edges**  
  - Reconstruction is "hole-free, but contains small bridges due to space carving"
  - Topological genus is larger than it appears
- **Size**: 543,652 vertices, 1,087,716 triangles (full resolution)
- **Format**: PLY files at multiple resolutions
- **Why This Matters**: This is a **widely cited benchmark** for mesh repair algorithms in academic research
- **Use Case**: Test non-manifold repair, self-intersection detection, topological analysis

### Dragon (`dragon_recon/`)
- **Source**: Stanford 3D Scanning Repository  
- **Known Defects**:
  - Contains **numerous small holes**
  - Incomplete surface from scanning artifacts
- **Size**: 566,098 vertices, 1,132,830 triangles
- **Format**: PLY files at multiple resolutions
- **Why This Matters**: Real-world hole detection test case from 3D scanning
- **Use Case**: Test hole detection, boundary loop extraction, hole filling

## Synthetic Defective Meshes (`defective/`)

Small, controlled test cases with specific defects for unit testing:

### `duplicate_vertices.obj`
- **Defect**: Exact duplicate vertices at same position
- **Why**: Tests basic vertex welding
- **Expected Fix**: Merge vertices 0 and 3

### `near_duplicate_vertices.obj`
- **Defect**: Vertices within tolerance distance (0.0001 units apart)
- **Why**: Tests tolerance-based vertex welding  
- **Expected Fix**: Merge with tolerance >= 0.001

### `t_junction.obj`
- **Defect**: T-junction where edges don't connect properly
- **Why**: Common in CAD exports, boolean operations
- **Expected Fix**: Close gaps by welding nearby vertices

### `non_manifold_edge.obj`
- **Defect**: Edge shared by 3 faces (violates 2-manifold property)
- **Why**: Tests non-manifold edge detection and repair
- **Expected Fix**: Split edge so each has ≤2 incident faces

### `non_manifold_vertex.obj`
- **Defect**: Bowtie configuration - vertex shared by disconnected fans
- **Why**: Tests non-manifold vertex detection
- **Expected Fix**: Split vertex for each disconnected fan

### `hole_mesh.obj`
- **Defect**: Missing faces create boundary edges
- **Why**: Tests hole detection
- **Expected Fix**: Fill hole with new faces

## Severely Broken Meshes (`mesh-repair-test-models/`)

From [caretdashcaret/MeshRepairTestModels](https://github.com/caretdashcaret/MeshRepairTestModels) - meshes that are so broken they can't even be previewed:

### `models/double_cube.stl`
- **Defect**: "has obvious holes that need to be filled"
- **Size**: ~1KB
- **Format**: STL
- **Purpose**: Simple hole filling test

### `models/rakdos.stl`
- **Defect**: "has many tiny non-manifold holes"  
- **Size**: ~10MB
- **Corruption Level**: **Severe** - GitHub 3D preview cannot render it
- **Format**: STL
- **Purpose**: Stress test for extreme mesh corruption
- **License**: CC BY-NC-SA 3.0 (Jenny/CaretDashCaret)

## Test Strategy

### 1. Demonstrate Defects EXIST (Before Repair)
**Goal**: Prove these meshes are actually broken

```rust
#[test]
fn test_happy_buddha_has_defects() {
    let mesh = load_ply("happy_recon/happy_vrip_res2.ply");
    
    // Should NOT be manifold
    assert_ne!(mesh.is_manifold(), RedgeManifoldness::IsManifold);
    
    // Count non-manifold edges
    let non_manifold_edges = mesh.meta_edges()
        .filter(|e| e.hedge().radial_loop().count() > 2)
        .count();
    
    assert!(non_manifold_edges > 0, "Happy Buddha should have non-manifold edges");
}
```

### 2. Show Repair WORKS (After Repair)

```rust
#[test]
fn test_vertex_welding_fixes_duplicates() {
    let mesh = load_obj("defective/duplicate_vertices.obj");
    
    // BEFORE: Has 4 vertices (including duplicate)
    assert_eq!(mesh.vert_count(), 4);
    
    let (repaired, report) = mesh.weld_vertices(1e-10);
    
    // AFTER: Has 3 vertices
    assert_eq!(repaired.vert_count(), 3);
    assert!(report.vertices_merged > 0);
}
```

### 3. Validate on Real-World Data

Use Happy Buddha and Dragon to prove algorithms work on actual scanned geometry with real defects, not just toy examples.

## Why These Specific Meshes?

### Academic Citations
- **Happy Buddha**: Used in mesh repair papers since it's a known benchmark
- **Dragon**: Classic test case for hole detection from Stanford
- **Rakdos**: Real-world severely corrupted mesh

### Documented Defects
We're not guessing - these defects are **explicitly documented**:
- Stanford repo states: "contains small bridges", "numerous small holes"  
- Happy Buddha paper: "153 self-intersections, 5581 non-manifold edges"
- CaretDashCaret repo: "has many tiny non-manifold holes"

### Real Problems
These aren't artificial - they're from:
- 3D scanning artifacts (Happy Buddha, Dragon)
- Space carving reconstruction issues (Happy Buddha)
- User-created 3D printing models (Rakdos)

## File Format Notes

- **PLY**: Stanford models (need PLY loader or conversion)
- **OBJ**: Synthetic defects (already supported via wavefront_loader)
- **STL**: MeshRepairTestModels (need STL loader or conversion)

## Converting PLY to OBJ

If needed, convert using:
```bash
# Using meshlabserver
meshlabserver -i happy_vrip_res2.ply -o happy_vrip_res2.obj

# Or Python with trimesh
python3 -c "import trimesh; trimesh.load('happy_vrip_res2.ply').export('happy_vrip_res2.obj')"
```

## References

1. **Thingi10K Dataset** - Zhou, Q. and Jacobson, A. (2016)
   - 10,000 models with real-world defects from Thingiverse
   - Available: https://github.com/Thingi10K/Thingi10K

2. **Stanford 3D Scanning Repository** - Curless, B. and Levoy, M.
   - Classic scanned models with documented issues
   - https://graphics.stanford.edu/data/3Dscanrep/

3. **Mesh repairing using topology graphs** - Li, X. et al. (2021)
   - Journal of Computational Design and Engineering
   - Uses Happy Buddha as benchmark (153 self-intersections, 5581 non-manifold edges)

## What We DON'T Include

❌ Clean, perfect meshes (Stanford Bunny's "clean" version, Utah Teapot, etc.)  
❌ Models without documented defects  
❌ Synthetic "nice" geometry  

## What We DO Include

✅ Meshes with **documented, cited defects**  
✅ Real-world broken geometry from scanning/reconstruction  
✅ Severely corrupted meshes that demonstrate necessity of repair  
✅ Specific, reproducible test cases for each defect type
