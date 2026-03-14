use std::env::args;
use std::path::PathBuf;
use colored::Colorize;

use anyhow::Result;
use log::info;
use mesh_splitting::mesh::Mesh;
use mesh_splitting::primitives::PrimitiveIdx;
use obj::Obj;

use mesh_splitting::mesh::parse_obj;
use mesh_splitting::Save;
use mesh_splitting::Collide;
use mesh_splitting::Split;


fn main() -> Result<()> {
    env_logger::init();

    let obj = Obj::load("./test/Remesh.obj")?;

    println!("{:?}", obj);

    let (mut meshes, verts, norms, faces) = parse_obj(obj);
    //let (meshes, ..) = parse_obj(obj);

    println!("Verts: {}", verts.borrow().len());
    println!("Norms: {}", norms.borrow().len());
    println!("Faces: {}", faces.borrow().len());
    for mesh in &meshes {
        println!(" > Mesh: {}", mesh.name.green());
        println!("Mesh: {:?}", mesh.polygons);
    }

    info!("Start splitting meshes");

    let mut resolved_meshes = Vec::new();
    let mut idx = meshes.len();
    for i in 0..meshes.len() {
        for j in i+1..meshes.len() {
            let mesh_a = &meshes[i];
            let mesh_b = &meshes[j];
            info!("Mesh {} and {} overlap: {}", mesh_a.name, mesh_b.name, mesh_a.overlap(mesh_b));
            let inter = mesh_a.intersect(mesh_b);
            info!("Mesh {} and {} intersection: {:?}", mesh_a.name, mesh_b.name, inter);
            let (new_mesh_a, inter) = mesh_a.split(inter);
            let (new_mesh_b, inter) = mesh_b.split(inter);
            //if !inter.is_empty() {
            //    let inter_mesh = Mesh {
            //        name: format!("{}-{}", mesh_a.name, mesh_b.name),
            //        idx: PrimitiveIdx::Global(idx),
            //        polygons: inter.into(),
            //        verts: mesh_a.verts.clone(),
            //        norms: mesh_a.norms.clone(),
            //        faces: mesh_a.faces.clone(),
            //    };
            //    resolved_meshes.push(inter_mesh);
            //    idx += 1;
            //}
            meshes[i] = new_mesh_a;
            meshes[j] = new_mesh_b;
        }
        resolved_meshes.push(meshes[i].clone());
    }

    let output_dir: PathBuf = args().nth(1).unwrap_or("./test/".to_string()).into();

    println!("Re-export the mesh to re-export_meshes.obj");
    resolved_meshes.save(&output_dir.join("re-export_meshes.obj"))?;

    Ok(())
}
