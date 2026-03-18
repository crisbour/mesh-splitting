use std::env::args;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use colored::Colorize;

use anyhow::Result;
use obj::Obj;

use mesh_splitting::mesh::parse_obj;
use mesh_splitting::mesh::remesh;
use mesh_splitting::Save;

fn main() -> Result<()> {
    env_logger::init();

    let obj_filepath = args().nth(1).unwrap_or("./test/WaterTank.obj".to_string());

    //let obj = Obj::load("./test/Remesh.obj")?;
    let obj = Obj::load(&obj_filepath)?;

    let now = Instant::now();

    let (meshes, verts, norms, faces) = parse_obj(obj.data);
    //let (meshes, ..) = parse_obj(obj);

    println!("Verts: {}", verts.borrow().len());
    println!("Norms: {}", norms.borrow().len());
    println!("Faces: {}", faces.borrow().len());
    for mesh in &meshes {
        println!(" > Mesh: {}", mesh.name.green());
        //println!("Mesh: {:?}", mesh.polygons);
    }

    let resolved_meshes = remesh(meshes)?;

    println!("Remeshed in {} seconds", now.elapsed().as_secs_f64());

   let output_dir: PathBuf = Path::new(&obj_filepath)
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));

    println!("Re-export the mesh to re-export_meshes.obj");
    resolved_meshes.save(&output_dir.join("re-export_meshes.obj"))?;

    Ok(())
}
