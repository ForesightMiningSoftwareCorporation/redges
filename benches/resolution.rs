use criterion::{criterion_group, criterion_main, Criterion};
use nalgebra::Vector3;
use quadric_simplification::StopCondition;
use redges::algorithms::remeshing::incremental_refinement_with_context;
use redges::algorithms::remeshing::RemeshingContext;
use redges::algorithms::remeshing::RemeshingParametersWithoutCollapse;
use redges::quadric_simplification::{QuadricSimplificationConfig, SimplificationStrategy};
use redges::*;
use spade::Triangulation;
use std::f64::consts::PI;

type Vec3 = Vector3<f64>;

pub fn mesh_from_polyline_and_elevations(
    mut polyline: Vec<Vec3>,
    top: f64,
    bot: f64,
) -> (Vec<Vec3>, Vec<u32>) {
    // Triangulate the bottom face.
    let edges: Vec<_> = (0..polyline.len())
        .map(|i| (i, (i + 1) % polyline.len()))
        .collect();

    let mut triangulation = spade::ConstrainedDelaunayTriangulation::<spade::Point2<_>>::default();

    let spade_points: Vec<_> = polyline
        .iter()
        .map(|p| triangulation.insert(spade::Point2::new(p.x, p.y)).unwrap())
        .collect();

    for (i, j) in edges {
        triangulation.add_constraint_and_split(spade_points[i], spade_points[j], |v| v);
    }

    let mut topology = Vec::<u32>::new();
    for face in triangulation.inner_faces() {
        let center_x = face.positions()[0].x + face.positions()[1].x + face.positions()[2].x;
        let center_y = face.positions()[0].y + face.positions()[1].y + face.positions()[2].y;

        let center = nalgebra::Point2::<f64>::new(center_x / 3., center_y / 3.);

        // Add to the mesh only those faces which are inside the polyline.
        let inside = point_in_poly2d(
            &center,
            &polyline
                .iter()
                .map(|p| nalgebra::Point2::<f64>::new(p.x, p.y))
                .collect::<Vec<_>>(),
        );
        if inside {
            topology.push(face.vertices()[0].index() as u32);
            topology.push(face.vertices()[1].index() as u32);
            topology.push(face.vertices()[2].index() as u32);
        }
    }

    for v in &mut polyline {
        v[2] = top;
    }

    // duplicate the top mesh and move it down.
    let n = polyline.len();
    for i in 0..n {
        let mut v = polyline[i];
        v[2] = bot;
        polyline.push(v);
    }

    // Add topology for the second face.
    let m = topology.len();
    for i in (0..m).rev() {
        let index = topology[i];
        topology.push(index + n as u32);
    }

    for i in 0..n {
        let top_left = i;
        let top_right = (i + 1) % n;
        let bot_left = n + i;
        let bot_right = n + ((i + 1) % n);

        // triangle 1.
        topology.push(bot_left as u32);
        topology.push(top_right as u32);
        topology.push(top_left as u32);

        // triangle 2.
        topology.push(bot_left as u32);
        topology.push(bot_right as u32);
        topology.push(top_right as u32);
    }

    (polyline, topology)
}

pub fn point_in_poly2d(pt: &nalgebra::Point2<f64>, poly: &[nalgebra::Point2<f64>]) -> bool {
    if poly.is_empty() {
        return false;
    }

    let mut winding = 0i32;

    for (i, a) in poly.iter().enumerate() {
        let b = poly[(i + 1) % poly.len()];
        let seg_dir = b - a;
        let dpt = pt - a;
        let perp = dpt.perp(&seg_dir);
        winding += match (dpt.y >= 0.0, b.y > pt.y) {
            (true, true) if perp < 0.0 => 1,
            (false, false) if perp > 0.0 => 1,
            _ => 0,
        };
    }

    winding % 2 == 1
}

pub fn subdivide_mesh(
    mesh: (Vec<Vec3>, Vec<usize>),
    number_of_subdivisions: usize,
    target_length: f64,
) -> (Vec<Vec3>, Vec<Vec<u32>>) {
    let redge =
        Redge::<(_, _, _)>::new(mesh.0, (), (), mesh.1.chunks(3).map(|c| c.iter().copied()));

    let mut context = RemeshingContext::from_mesh(redge);

    incremental_refinement_with_context(
        &mut context,
        RemeshingParametersWithoutCollapse {
            target_additional_vertices: number_of_subdivisions,
            ..Default::default()
        },
        |_i, _j| target_length,
        |&v| v,
    );

    let (vs, fs, _) = context.mesh.to_face_list();
    let fs = fs
        .into_iter()
        .map(|l| l.into_iter().map(|i| i as u32).collect())
        .collect();
    (vs, fs)
}

fn subdivision_benchmark(c: &mut Criterion) {
    let res = 50;
    let parametric_bean: Vec<_> = (0..res)
        .map(|i| {
            let t = (i as f64 / res as f64) * 2. * PI - PI;
            Vec3::new(t.cos() - (1. / (1. + 4. * t * t)), t.sin(), 0.)
        })
        .collect();

    let (vs, fs) = mesh_from_polyline_and_elevations(parametric_bean, 1.0, -1.0);

    let fs: Vec<_> = fs.iter().map(|i| *i as usize).collect();

    c.bench_function("redge subdivision", |b| {
        b.iter(|| {
            let _ = subdivide_mesh((vs.clone(), fs.clone()), 1_000_000, 0.4);
        });
    });
}

fn custom_criterion() -> Criterion {
    Criterion::default()
        .measurement_time(std::time::Duration::from_secs(300))
        .warm_up_time(std::time::Duration::from_secs(10))
}

fn simplification_benchmark(c: &mut Criterion) {
    let res = 1000;
    let parametric_bean: Vec<_> = (0..res)
        .map(|i| {
            let t = (i as f64 / res as f64) * 2. * PI - PI;
            Vec3::new(t.cos() - (1. / (1. + 4. * t * t)), t.sin(), 0.)
        })
        .collect();

    let (vs, fs) = mesh_from_polyline_and_elevations(parametric_bean, 1.0, -1.0);

    let fs: Vec<_> = fs.iter().map(|i| *i as usize).collect();

    c.bench_function("redge benchmark", |b| {
        b.iter(|| {
            let redge = Redge::<(_, _, _)>::new(
                vs.clone(),
                (),
                (),
                fs.chunks(3).map(|f| f.iter().copied()),
            );

            let face_count_before = redge.face_count();
            let redge = redge.clean_overlapping_faces();
            let (_redge, _) = quadric_simplification::quadric_simplify(
                redge,
                QuadricSimplificationConfig {
                    strategy: SimplificationStrategy::Conservative,
                    attribute_simplification:
                        quadric_simplification::AttributeSimplification::NoAttributeSimplification,
                    stop_condition: StopCondition::FaceCount(face_count_before / 10),
                },
                |_, _| false,
            );
        });
    });
}

criterion_group! {
    name = benches;
    config = custom_criterion();
    targets = subdivision_benchmark, simplification_benchmark
}
criterion_main!(benches);
