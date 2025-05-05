use redges::container_trait::FaceAttributeGetter;
use redges::quadric_simplification::{QuadricSimplificationConfig, SimplificationStrategy};
use redges::VertId;
use redges::{quadric_simplification, Redge};

use redges::wavefront_loader::*;
use std::fs::File;
use std::io::Write;

pub type Vec2 = nalgebra::Vector2<f64>;
pub type Vec3 = nalgebra::Vector3<f64>;

#[derive(Debug, Default, Clone)]
pub struct FaceData {
    pub uvs: [Vec2; 3],
    pub verts: [VertId; 3],
}

impl FaceAttributeGetter<f64> for FaceData {
    fn attribute(&self, vert_index: usize, attribute_id: usize) -> f64 {
        self.uvs[vert_index][attribute_id]
    }

    fn attribute_count(&self) -> usize {
        2
    }

    fn attribute_mut(&mut self, vert_index: usize, attribute_id: usize) -> &mut f64 {
        &mut self.uvs[vert_index][attribute_id]
    }

    fn inner_index(&self, vid: VertId) -> usize {
        self.verts
            .iter()
            .position(|id| *id == vid)
            .unwrap_or_else(|| {
                panic!(
                    "Asked for position of vertex {:?} in face {:?}",
                    vid, self.verts
                )
            })
    }

    fn attribute_vertices(&self) -> &[VertId] {
        &self.verts
    }

    fn attribute_vertices_mut(&mut self) -> &mut [VertId] {
        &mut self.verts
    }
}

fn export_to_obj(vertices: &[Vec3], faces: &Vec<FaceData>, path: &str) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for v in vertices {
        writeln!(file, "v {} {} {}", v.x, v.y, v.z,)?;
    }

    for face in faces {
        writeln!(file, "vt {} {}", face.uvs[0].x, face.uvs[0].y)?;
        writeln!(file, "vt {} {}", face.uvs[1].x, face.uvs[1].y)?;
        writeln!(file, "vt {} {}", face.uvs[2].x, face.uvs[2].y)?;
    }

    for (i, face) in faces.iter().enumerate() {
        writeln!(
            file,
            "f {}/{} {}/{} {}/{}",
            face.verts[0].to_index() + 1,
            i * 3 + 1, // vt
            face.verts[1].to_index() + 1,
            i * 3 + 2, // vt
            face.verts[2].to_index() + 1,
            i * 3 + 3, // vt
        )?;
    }

    Ok(())
}

fn export_to_obj_indices(
    vertices: &[Vec3],
    faces: &[Vec<usize>],
    path: &str,
) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for v in vertices {
        writeln!(file, "v {} {} {}", v.x, v.y, v.z,)?;
    }

    for face in faces.iter() {
        writeln!(file, "f {} {} {}", face[0] + 1, face[1] + 1, face[2] + 1,)?;
    }

    Ok(())
}

fn main() {
    let mut obj_data = ObjData::from_disk_file("assets/armadillo.obj");
    if obj_data.uv_face_indices.is_empty() {
        obj_data.uvs = vec![nalgebra::Vector2::new(0., 0.)];
        obj_data.uv_face_indices = vec![vec![0; 3]; obj_data.vertex_face_indices.len()];
    }

    assert!(obj_data.uv_face_indices.len() == obj_data.vertex_face_indices.len());
    let mut faces: Vec<_> = (0..obj_data.vertex_face_indices.len())
        .map(|i| FaceData {
            uvs: [
                obj_data.uvs[obj_data.uv_face_indices[i][0] as usize]
                    .clone()
                    .map(|s| s as f64),
                obj_data.uvs[obj_data.uv_face_indices[i][1] as usize]
                    .clone()
                    .map(|s| s as f64),
                obj_data.uvs[obj_data.uv_face_indices[i][2] as usize]
                    .clone()
                    .map(|s| s as f64),
            ],
            verts: [
                VertId(obj_data.vertex_face_indices[i][0] as usize),
                VertId(obj_data.vertex_face_indices[i][1] as usize),
                VertId(obj_data.vertex_face_indices[i][2] as usize),
            ],
        })
        .collect();

    if faces.is_empty() {
        for i in 0..obj_data.uv_face_indices.len() {
            faces.push(FaceData {
                uvs: [nalgebra::Vector2::default(); 3],
                verts: [
                    VertId(obj_data.uv_face_indices[i][0] as usize),
                    VertId(obj_data.uv_face_indices[i][1] as usize),
                    VertId(obj_data.uv_face_indices[i][2] as usize),
                ],
            });
        }
    }

    let vertices: Vec<_> = obj_data
        .vertices
        .iter()
        .map(|v| v.map(|s| s as f64))
        .collect();

    let indices: Vec<_> = obj_data
        .vertex_face_indices
        .iter()
        .map(|l| {
            l.clone()
                .into_iter()
                .map(|i| i as usize)
                .collect::<Vec<_>>()
        })
        .collect();

    for _ in 0..1000 {
        let (vs, fs) = simplify_example(vertices.clone(), faces.clone(), indices.clone());
    }
    // export_to_obj(&vs, &fs, "tmp/simple_with_attributes.obj").unwrap();
    //
    // let (vs, indices) = simplify_example_geometry_only(vertices.clone(), indices.clone());
    // export_to_obj_indices(&vs, &indices, "tmp/simple_without_attributes.obj").unwrap();
}

fn simplify_example(
    vertices: Vec<Vec3>,
    faces: Vec<FaceData>,
    indices: Vec<Vec<usize>>,
) -> (Vec<Vec3>, Vec<FaceData>) {
    let redge = Redge::<(_, _, _)>::new(
        vertices,
        (),
        faces,
        indices.iter().map(|f| f.iter().copied()),
    );

    let face_count_before = redge.face_count();
    let redge = redge.clean_overlapping_faces();
    let (redge, _) = quadric_simplification::quadric_simplify(
        redge,
        QuadricSimplificationConfig {
            strategy: SimplificationStrategy::Conservative,
            attribute_simplification:
                quadric_simplification::AttributeSimplification::SimplifyAtributes,
            target_face_count: face_count_before / 10,
        },
        |_, _| false,
    );

    let (vs, _ids, fs) = redge.to_face_list();

    (vs, fs)
}

fn simplify_example_geometry_only(
    vertices: Vec<Vec3>,
    indices: Vec<Vec<usize>>,
) -> (Vec<Vec3>, Vec<Vec<usize>>) {
    let redge =
        Redge::<(_, _, _)>::new(vertices, (), (), indices.iter().map(|f| f.iter().copied()));

    let face_count_before = redge.face_count();

    let redge = redge.clean_overlapping_faces();
    let (redge, _) = quadric_simplification::quadric_simplify(
        redge,
        QuadricSimplificationConfig {
            strategy: SimplificationStrategy::Conservative,
            attribute_simplification:
                quadric_simplification::AttributeSimplification::NoAttributeSimplification,
            target_face_count: face_count_before / 10,
        },
        |_, _| false,
    );

    let (vs, ids, _fs) = redge.to_face_list();

    (vs, ids)
}
