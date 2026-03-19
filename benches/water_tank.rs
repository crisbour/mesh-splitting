use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use obj::Obj;
use mesh_splitting::mesh::parse_obj;
use mesh_splitting::mesh::remesh;

fn criterion_benchmark(c: &mut Criterion) {
    // Benchmark loading
    c.bench_function("load_obj", |b| {
        b.iter(|| {
            let obj = Obj::load(black_box("./test/WaterTank.obj")).expect("Failed to load OBJ file");
            black_box(obj);
        });
    });

    // Benchmark parsing
    let obj = Obj::load("./test/WaterTank.obj").expect("Failed to load OBJ file");
    c.bench_function("parse_obj", |b| {
        b.iter(|| {
            let (meshes, verts, norms, faces) = parse_obj(black_box(&obj.data));
            black_box((meshes, verts, norms, faces));
        });
    });

    // Benchmark remeshing
    let (meshes, verts, norms, faces) = parse_obj(&obj.data);
    c.bench_function("remesh", |b| {
        b.iter(|| {
            let result = remesh(black_box(meshes.clone()));
            assert!(result.is_ok());
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
