//! Helper functions to load various mesh formats for testing
//!
//! This module provides loaders for:
//! - OBJ (via existing wavefront_loader)
//! - PLY (via ply-rs)
//! - STL (via stl_io)

use nalgebra::Vector3;
use ply_rs::{parser, ply};
use redges::{wavefront_loader::ObjData, Redge};
use std::fs::File;
use std::io::BufReader;

type TestMesh = Redge<(Vec<Vector3<f64>>, (), ())>;

/// Load OBJ file
pub fn load_obj_mesh(path: &str) -> TestMesh {
    let obj_data = ObjData::from_disk_file(path);
    let vertices: Vec<_> = obj_data
        .vertices
        .iter()
        .map(|v| v.map(|s| s as f64))
        .collect();
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

/// Load PLY file
pub fn load_ply_mesh(path: &str) -> TestMesh {
    let file = File::open(path).expect("Failed to open PLY file");
    let mut buf_reader = BufReader::new(file);

    // Create parser
    let parser = parser::Parser::<ply::DefaultElement>::new();

    // Read the PLY file
    let ply = parser
        .read_ply(&mut buf_reader)
        .expect("Failed to parse PLY file");

    // Extract vertices
    let mut vertices = Vec::new();
    if let Some(vertex_list) = ply.payload.get("vertex") {
        for vertex_elem in vertex_list {
            let x = match vertex_elem.get("x") {
                Some(ply::Property::Float(v)) => *v as f64,
                Some(ply::Property::Double(v)) => *v,
                _ => 0.0,
            };
            let y = match vertex_elem.get("y") {
                Some(ply::Property::Float(v)) => *v as f64,
                Some(ply::Property::Double(v)) => *v,
                _ => 0.0,
            };
            let z = match vertex_elem.get("z") {
                Some(ply::Property::Float(v)) => *v as f64,
                Some(ply::Property::Double(v)) => *v,
                _ => 0.0,
            };
            vertices.push(Vector3::new(x, y, z));
        }
    }

    // Extract faces
    let mut faces = Vec::new();
    if let Some(face_list) = ply.payload.get("face") {
        for face_elem in face_list {
            if let Some(ply::Property::ListInt(indices)) = face_elem.get("vertex_indices") {
                let face_indices: Vec<usize> = indices.iter().map(|&i| i as usize).collect();
                faces.push(face_indices);
            } else if let Some(ply::Property::ListUInt(indices)) = face_elem.get("vertex_indices")
            {
                let face_indices: Vec<usize> = indices.iter().map(|&i| i as usize).collect();
                faces.push(face_indices);
            }
        }
    }

    Redge::new(
        vertices,
        (),
        (),
        faces.iter().map(|f| f.iter().copied()),
    )
}

/// Load STL file
pub fn load_stl_mesh(path: &str) -> TestMesh {
    let mut file = File::open(path).expect("Failed to open STL file");
    let stl = stl_io::read_stl(&mut file).expect("Failed to parse STL file");

    // stl_io returns an IndexedMesh with vertices and faces
    let vertices: Vec<_> = stl
        .vertices
        .iter()
        .map(|v| Vector3::new(v[0] as f64, v[1] as f64, v[2] as f64))
        .collect();

    let faces: Vec<_> = stl
        .faces
        .iter()
        .map(|f| vec![f.vertices[0], f.vertices[1], f.vertices[2]])
        .collect();

    Redge::new(
        vertices,
        (),
        (),
        faces.iter().map(|f| f.iter().copied()),
    )
}
