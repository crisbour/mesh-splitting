use gungraun::{main, library_benchmark_group, library_benchmark, LibraryBenchmarkConfig};
use std::hint::black_box;

use obj::Obj;
use mesh_splitting::mesh::parse_obj;
use mesh_splitting::mesh::remesh;

#[library_benchmark]
#[benches::multiple("test/WaterTank.obj", "test/Remesh.obj")]
fn bench_remesh_water_tank(path: &str) {
    let path = format!("{}/{}", env!("CARGO_MANIFEST_DIR"), path);
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
