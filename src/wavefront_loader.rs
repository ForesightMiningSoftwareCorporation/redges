use std::fmt::{Debug, Display};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::ops::{AddAssign, Index};

use linear_isomorphic::prelude::*;

type Vec3 = [f32; 3];
type Vec2 = [f32; 2];

#[derive(Debug)]
pub struct ObjData {
    pub vertices: Vec<Vec3>,
    pub normals: Vec<Vec3>,
    pub uvs: Vec<Vec2>,
    pub vertex_face_indices: Vec<Vec<u64>>,
    pub normal_face_indices: Vec<Vec<u64>>,
    pub uv_face_indices: Vec<Vec<u64>>,
    pub lines: Vec<[u64; 2]>,
    pub objects: Vec<ObjectRanges>,
}

#[derive(Debug)]
pub struct ObjectRanges {
    pub name: String,
    pub positions: (u64, u64),
    pub normals: (u64, u64),
    pub uvs: (u64, u64),
    pub lines: (u64, u64),
    pub polygons: [(u64, u64); 3],
}

impl ObjData {
    pub fn from_disk_file(path: &str) -> Self {
        let file = File::open(path).expect(&format!("Cannot find {}.", path));
        let reader = BufReader::new(file);

        let mut vert_topology = Vec::<Vec<u64>>::new();
        let mut normal_topology = Vec::<Vec<u64>>::new();
        let mut uv_topology = Vec::<Vec<u64>>::new();
        let mut lines = Vec::<[u64; 2]>::new();

        let mut vertices = Vec::<Vec3>::new();
        let mut normals = Vec::<Vec3>::new();
        let mut uvs = Vec::<Vec2>::new();
        let mut object_ranges = Vec::<ObjectRanges>::new();
        for line in reader.lines() {
            let line_str = line.unwrap();
            let line_str = line_str.trim();
            let tokens = line_str.split(" ").collect::<Vec<&str>>();
            match tokens[0] {
                "o" => {
                    object_ranges.push(ObjectRanges {
                        name: tokens[1..].join(" "),
                        positions: (vertices.len() as u64, 0),
                        normals: (normals.len() as u64, 0),
                        uvs: (uvs.len() as u64, 0),
                        lines: (lines.len() as u64, 0),
                        polygons: [
                            (vert_topology.len() as u64, 0),
                            (normal_topology.len() as u64, 0),
                            (uv_topology.len() as u64, 0),
                        ],
                    });
                }
                "v" => {
                    add_vec_3d(&tokens[1..], &mut vertices);
                }
                "vn" => {
                    add_vec_3d(&tokens[1..], &mut normals);
                }
                "vt" => {
                    add_vec_2d(&tokens[1..], &mut uvs);
                }
                "l" => {
                    add_line(&tokens[1..], &mut lines);
                }
                "f" => add_face(
                    &tokens[1..],
                    &mut vert_topology,
                    &mut normal_topology,
                    &mut uv_topology,
                ),
                _ => continue,
            }

            if let Some(or) = object_ranges.last_mut() {
                or.positions.1 = vertices.len() as u64;
                or.normals.1 = normals.len() as u64;
                or.uvs.1 = uvs.len() as u64;
                or.lines.1 = lines.len() as u64;

                or.polygons[0].1 = vert_topology.len() as u64;
                or.polygons[1].1 = normal_topology.len() as u64;
                or.polygons[2].1 = uv_topology.len() as u64;
            }
        }

        ObjData {
            vertices,
            normals,
            uvs,
            vertex_face_indices: vert_topology,
            normal_face_indices: normal_topology,
            uv_face_indices: uv_topology,
            lines,
            objects: object_ranges,
        }
    }

    pub fn export<'a, T>(mesh: &'a T, path: &str)
    where
        T: WaveFrontCompatible<'a>,
    {
        std::fs::create_dir_all(std::path::Path::new(path).parent().unwrap()).unwrap();
        let mut file = File::create(path).unwrap();

        let one = 1;
        for point in mesh.pos_iterator() {
            writeln!(file, "v {} {} {}", point[0], point[1], point[2]).unwrap();
        }

        for norm in mesh.norm_iterator() {
            writeln!(file, "vn {} {} {}", norm[0], norm[1], norm[2]).unwrap();
        }

        for uv in mesh.uv_iterator() {
            writeln!(file, "vt {} {}", uv[0], uv[1]).unwrap();
        }

        for line in mesh.segment_iterator() {
            writeln!(file, "l {} {}", line[0] + one, line[1] + one).unwrap();
        }

        let mut indices = Vec::new();
        let mut vert_face_count = 0;
        for face in mesh.pos_index_iterator() {
            vert_face_count += 1;

            indices.push(Vec::new());
            let face_ids: &mut Vec<_> = indices.last_mut().unwrap();

            for pos_id in face {
                face_ids.push([Some(pos_id), None, None]);
            }
        }

        let mut uv_face_count = 0;
        for face in mesh.uv_index_iterator() {
            for (local_id, uv_id) in face.enumerate() {
                indices[uv_face_count][local_id][1] = Some(uv_id);
            }
            uv_face_count += 1;
        }

        debug_assert!(uv_face_count == 0 || uv_face_count == vert_face_count);

        let mut norm_face_count = 0;
        for face in mesh.norm_index_iterator() {
            for (local_id, norm_id) in face.enumerate() {
                indices[norm_face_count][local_id][2] = Some(norm_id);
            }
            norm_face_count += 1;
        }

        debug_assert!(norm_face_count == 0 || norm_face_count == vert_face_count);

        for face in indices {
            let mut face_data = "f ".to_string();

            for [pos_id, uv_id, norm_id] in face {
                let pos_str = (pos_id.unwrap() + one).to_string();
                let uv_str = match uv_id {
                    None => "".to_string(),
                    Some(x) => (x + one).to_string(),
                };
                let norm_str = match norm_id {
                    None => "".to_string(),
                    Some(x) => (x + one).to_string(),
                };

                face_data.push_str(format!("{}/{}/{} ", pos_str, uv_str, norm_str).as_str());
            }

            face_data.push('\n');

            write!(file, "{}", face_data).unwrap();
        }
    }
}

pub trait WaveFrontCompatible<'a> {
    type Scalar: num_traits::Float + Debug + Display;

    fn pos_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]>;
    fn uv_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 2]>;
    fn norm_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]>;

    fn segment_iterator(&'a self) -> impl Iterator<Item = [usize; 2]>;
    fn pos_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>>;
    fn uv_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>>;
    fn norm_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>>;
}

impl<'a> WaveFrontCompatible<'a> for ObjData {
    type Scalar = f32;

    fn pos_iterator(&'a self) -> impl Iterator<Item = [f32; 3]> {
        self.vertices.iter().map(|v| [v[0], v[1], v[2]])
    }

    fn uv_iterator(&'a self) -> impl Iterator<Item = [f32; 2]> {
        self.uvs.iter().map(|v| [v[0], v[1]])
    }

    fn norm_iterator(&'a self) -> impl Iterator<Item = [f32; 3]> {
        self.normals.iter().map(|v| [v[0], v[1], v[2]])
    }

    fn segment_iterator(&'a self) -> impl Iterator<Item = [usize; 2]> {
        self.lines.iter().map(|[i, j]| [*i as usize, *j as usize])
    }

    fn pos_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        self.vertex_face_indices
            .iter()
            .map(|l| l.iter().map(|i| *i as usize))
    }

    fn uv_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        self.uv_face_indices
            .iter()
            .map(|l| l.iter().map(|i| *i as usize))
    }

    fn norm_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        self.normal_face_indices
            .iter()
            .map(|l| l.iter().map(|i| *i as usize))
    }
}

impl<'a, V, S> WaveFrontCompatible<'a> for Vec<V>
where
    V: Index<usize, Output = S>,
    S: num_traits::Float + Display + Debug,
{
    type Scalar = S;

    fn pos_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]> {
        self.iter().map(|v| [v[0], v[1], v[2]])
    }

    fn uv_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 2]> {
        std::iter::empty()
    }

    fn norm_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]> {
        std::iter::empty()
    }

    fn segment_iterator(&'a self) -> impl Iterator<Item = [usize; 2]> {
        std::iter::empty()
    }

    fn pos_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }

    fn uv_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }

    fn norm_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }
}

impl<'a, V, S> WaveFrontCompatible<'a> for (&Vec<V>, &Vec<usize>)
where
    V: Index<usize, Output = S>,
    S: num_traits::Float + Display + Debug,
{
    type Scalar = S;

    fn pos_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]> {
        self.0.iter().map(|v| [v[0], v[1], v[2]])
    }

    fn uv_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 2]> {
        std::iter::empty()
    }

    fn norm_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]> {
        std::iter::empty()
    }

    fn segment_iterator(&'a self) -> impl Iterator<Item = [usize; 2]> {
        std::iter::empty()
    }

    fn pos_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        debug_assert!(self.1.len() % 3 == 0);
        self.1.chunks(3).map(|chunk| chunk.iter().copied())
    }

    fn uv_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }

    fn norm_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }
}

impl<'a, V, S, I> WaveFrontCompatible<'a> for (&Vec<V>, &Vec<Vec<I>>)
where
    V: Index<usize, Output = S>,
    S: num_traits::Float + Display + Debug,
    usize: TryFrom<I>,
    I: num_traits::PrimInt + std::fmt::Display,
{
    type Scalar = S;

    fn pos_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]> {
        self.0.iter().map(|v| [v[0], v[1], v[2]])
    }

    fn uv_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 2]> {
        std::iter::empty()
    }

    fn norm_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]> {
        std::iter::empty()
    }

    fn segment_iterator(&'a self) -> impl Iterator<Item = [usize; 2]> {
        std::iter::empty()
    }

    fn pos_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        self.1.iter().map(|face| {
            face.iter()
                .map(|i| usize::try_from(*i).unwrap_or(usize::MAX))
        })
    }

    fn uv_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }

    fn norm_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }
}

impl<'a, V, S, I> WaveFrontCompatible<'a> for (&Vec<V>, &Vec<[I; 2]>)
where
    V: VectorSpace<Scalar = S>,
    S: num_traits::Float + Display + Debug,
    usize: TryFrom<I>,
    I: num_traits::PrimInt + Display,
{
    type Scalar = S;

    fn pos_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]> {
        self.0.iter().map(|v| [v[0], v[1], v[2]])
    }

    fn uv_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 2]> {
        std::iter::empty()
    }

    fn norm_iterator(&'a self) -> impl Iterator<Item = [Self::Scalar; 3]> {
        std::iter::empty()
    }

    fn segment_iterator(&'a self) -> impl Iterator<Item = [usize; 2]> {
        self.1.iter().map(|[i, j]| {
            [
                usize::try_from(*i).unwrap_or(usize::MAX),
                usize::try_from(*j).unwrap_or(usize::MAX),
            ]
        })
    }

    fn pos_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }

    fn uv_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }

    fn norm_index_iterator(&'a self) -> impl Iterator<Item = impl Iterator<Item = usize>> {
        let empty_iterator: std::iter::Empty<std::iter::Empty<usize>> = std::iter::empty();
        empty_iterator
    }
}

fn add_vec_3d(tokens: &[&str], vecs: &mut Vec<Vec3>) {
    let vec = [
        tokens[0].parse::<f32>().unwrap(),
        tokens[1].parse::<f32>().unwrap(),
        tokens[2].parse::<f32>().unwrap(),
    ];
    vecs.push(vec);
}

fn add_vec_2d(tokens: &[&str], vecs: &mut Vec<Vec2>) {
    let vec = [
        tokens[0].parse::<f32>().unwrap(),
        tokens[1].parse::<f32>().unwrap(),
    ];
    vecs.push(vec);
}

fn add_line(tokens: &[&str], vecs: &mut Vec<[u64; 2]>) {
    assert!(tokens.len() == 2, "have {} expected 2", tokens.len());
    let index_1 = tokens[0].parse::<u64>().unwrap() - 1;
    let index_2 = tokens[1].parse::<u64>().unwrap() - 1;

    vecs.push([index_1, index_2])
}

fn add_face(
    tokens: &[&str],
    vert_topology: &mut Vec<Vec<u64>>,
    normal_topology: &mut Vec<Vec<u64>>,
    uv_topology: &mut Vec<Vec<u64>>,
) {
    assert!(
        tokens.len() >= 3,
        "have {} out of a minimum of 3",
        tokens.len()
    );
    vert_topology.push(Vec::new());
    uv_topology.push(Vec::new());
    normal_topology.push(Vec::new());

    for token in tokens {
        let inner_tokens = token.split("/").collect::<Vec<&str>>();

        assert!(inner_tokens.len() <= 3 && inner_tokens.len() >= 1);

        vert_topology
            .last_mut()
            .unwrap()
            .push(inner_tokens[0].parse::<u64>().unwrap() - 1);

        // Only positions are specified.
        if inner_tokens.len() == 1 {
            continue;
        }

        // Only positions and texture coordinates are specified.
        if inner_tokens.len() == 2 {
            uv_topology
                .last_mut()
                .unwrap()
                .push(inner_tokens[1].parse::<u64>().unwrap() - 1);
            continue;
        }

        if !inner_tokens[1].is_empty() {
            uv_topology
                .last_mut()
                .unwrap()
                .push(inner_tokens[1].parse::<u64>().unwrap() - 1);
        }

        if !inner_tokens[2].is_empty() {
            normal_topology
                .last_mut()
                .unwrap()
                .push(inner_tokens[2].parse::<u64>().unwrap() - 1);
        }
    }

    if uv_topology.last().unwrap().is_empty() {
        uv_topology.pop();
    }
    if normal_topology.last().unwrap().is_empty() {
        normal_topology.pop();
    }
}
