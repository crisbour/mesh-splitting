use gungraun::{main, library_benchmark_group, library_benchmark, LibraryBenchmarkConfig};
use std::hint::black_box;

use obj::Obj;
use mesh_splitting::mesh::parse_obj;
use mesh_splitting::mesh::remesh;

#[library_benchmark]
fn bench_remesh_water_tank() {
    let path = format!("{}/test/WaterTank.obj", env!("CARGO_MANIFEST_DIR"));
    let obj = Obj::load(black_box(&path)).expect("Failed to load OBJ file");
    //let obj = Obj::load(black_box("./test/WaterTank.obj")).expect("Failed to load OBJ file");
    let (meshes, _verts, _norms, _faces) = parse_obj(obj.data);
    let result = remesh(black_box(meshes.clone()));
    assert!(result.is_ok());
}

library_benchmark_group!(
    name = remesh_group,
    benchmarks = [ bench_remesh_water_tank ]
);

main!(library_benchmark_groups = remesh_group);
//main!(remesh_group);
