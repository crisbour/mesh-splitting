use std::{cell::RefCell, rc::Rc};

use anyhow::Result;
use obj::Obj;

use mesh_splitting::mesh::parse_obj;


fn main() -> Result<()> {
    env_logger::init();

    let obj = Obj::load("./test/Remesh.obj")?;

    println!("{:?}", obj);

    let verts = Rc::new(RefCell::new(Vec::new()));
    let norms = Rc::new(RefCell::new(Vec::new()));
    let faces = Rc::new(RefCell::new(Vec::new()));

    let meshes = parse_obj(obj, verts.clone(), norms.clone(), faces.clone());

    println!("Verts: {:?}", verts.borrow());
    println!("Norms: {:?}", norms.borrow());
    println!("Faces: {:?}", faces.borrow());
    for mesh in &meshes {
        println!("Mesh: {:?}", mesh.polygons);
    }

    Ok(())
}
