use anyhow::Result;

use obj::{Group, IndexTuple, Obj, ObjData, Object, SimplePolygon};

use crate::{Save, mesh::Mesh, primitives::PrimitiveIdx};

fn vec_to_hash(vec: &[f32;3]) -> u128 {
    (vec[0].to_bits() as u128) << (0*32) |
    (vec[1].to_bits() as u128) << (1*32) |
    (vec[2].to_bits() as u128) << (2*32)
}

fn mesh_to_obj(mesh: &Mesh) -> Result<Object> {

    let mut object = Object::new(mesh.name.clone());
    let mut group = Group::new("default".to_string());
    let faces_ref = mesh.faces.borrow();

    group.polys = mesh
        .polygons
        .iter()
        .map(|tri_idx| match tri_idx {
            PrimitiveIdx::Local(_) => {
                panic!("Expected global indices for export, found: {:?}", tri_idx);
            },
            PrimitiveIdx::Global(idx) => *idx,
        })
        .map(|idx| &faces_ref[idx])
        .map(|tri| {
            let verts = tri.verts
                .map(|vert| match vert.idx {
                    PrimitiveIdx::Local(_) => {
                        panic!("Expected global indices for export, found: {:?}", tri.verts);
                    },
                    PrimitiveIdx::Global(idx) => idx,
                });
            let norms = tri.norms
                .map(|norms|
                    norms.map(|norm| match norm.idx {
                        PrimitiveIdx::Local(_) => {
                            panic!("Expected global indices for export, found: {:?}", norm.idx);
                        },
                        PrimitiveIdx::Global(idx) => idx,
                    })
                );
            let norm_options: [Option<usize>; 3] = match norms {
                Some(arr) => [Some(arr[0]), Some(arr[1]), Some(arr[2])],
                None => [None, None, None],
            };
            (verts, norm_options)
        })
        .map(|(verts, norms)| SimplePolygon(
            verts.iter().zip(norms).map(|(vert, norm)| IndexTuple(*vert, None, norm)).collect()
        ))
        .collect();

    object.groups.push(group);

    Ok(object)
}

impl Save for Vec<Mesh> {
    fn save_data(&self, path: &std::path::Path) -> Result<()> {
        let mut obj_data = ObjData::default();

        obj_data.position = self[0].verts
            .borrow()
            .iter()
            .map(|&vert| [vert.x, vert.y, vert.z].map(|p| p as f32))
            .collect();

        obj_data.normal = self[0].norms
            .borrow()
            .iter()
            .map(|&dir| [dir.x, dir.y, dir.z].map(|p| p as f32))
            .collect();

        for mesh in self {
            let object = mesh_to_obj(mesh)?;
            obj_data.objects.push(object);
        }

        let obj = Obj { data: obj_data, path: path.to_path_buf() };

        obj.save(path)?;
        Ok(())
    }
}

