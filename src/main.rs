use std::env::args;
use std::path::PathBuf;
use colored::Colorize;

use anyhow::Result;
use log::info;
use mesh_splitting::mesh::Mesh;
use mesh_splitting::mesh::remesh;
use mesh_splitting::primitives::PrimitiveIdx;
use obj::Obj;

use mesh_splitting::mesh::parse_obj;
use mesh_splitting::Save;
use mesh_splitting::Collide;
use mesh_splitting::Split;

fn main() -> Result<()> {
    env_logger::init();

    let obj = Obj::load("./test/Remesh.obj")?;

    //println!("{:?}", obj);

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

    let output_dir: PathBuf = args().nth(1).unwrap_or("./test/".to_string()).into();

    println!("Re-export the mesh to re-export_meshes.obj");
    resolved_meshes.save(&output_dir.join("re-export_meshes.obj"))?;

    Ok(())
}
