use std::{cell::RefCell, rc::Rc};

use nalgebra::{Point3, Unit, Vector3};
use obj::Obj;

use crate::{
    Collide, IdxTriangle, Split, SplitEdges,
    primitives::{IdxEdge, Normal, PrimitiveIdx, Vertex},
};

type Dir3 = Unit<Vector3<f64>>;

pub fn parse_obj(
    obj: Obj,
    verts: Rc<RefCell<Vec<Point3<f64>>>>,
    norms: Rc<RefCell<Vec<Dir3>>>,
    faces: Rc<RefCell<Vec<IdxTriangle>>>,
) -> Vec<Mesh> {
    *verts.borrow_mut() = obj
        .data
        .position
        .iter()
        .map(|vs| {
            let vs: Vec<f64> = vs.iter().map(|v| *v as f64).collect();
            Point3::new(vs[0], vs[1], vs[2])
        })
        .collect::<Vec<_>>();

    *norms.borrow_mut() = obj
        .data
        .normal
        .iter()
        .map(|vs| {
            let vs: Vec<f64> = vs.iter().map(|v| *v as f64).collect();
            Dir3::new_normalize(Vector3::new(vs[0], vs[1], vs[2]))
        })
        .collect::<Vec<_>>();

    obj.data
        .objects
        .iter()
        .enumerate()
        .map(|(idx, obj)| {
            // NOTE: Don't support multiple groups per object
            let obj_faces: Vec<IdxTriangle> = obj
                .groups
                .iter()
                .flat_map(|group| {
                    group.polys.iter().map(|poly| {
                        let tri_verts: Vec<_> = poly
                            .0
                            .iter()
                            .map(|idx_tuple| idx_tuple.0)
                            .map(|idx| Vertex::new(idx, &verts.borrow()))
                            .collect();
                        //let t_idx = poly.0.map(|idx_tuple| idx_tuple.1);
                        let tri_norms: Option<Vec<_>> = poly
                            .0
                            .iter()
                            .map(|idx_tuple| idx_tuple.2)
                            .map(|idx| idx.map(|i| Normal::new(i, &norms.borrow())))
                            .collect();
                        IdxTriangle::new(tri_verts, tri_norms)
                    })
                })
                .collect();

            let faces_len = faces.borrow().len();
            let polygons = (faces_len..(faces_len + obj_faces.len())).collect();
            faces.borrow_mut().extend(obj_faces);
            Mesh {
                name: obj.name.clone(),
                idx: PrimitiveIdx::Global(idx),
                polygons,
                verts: verts.clone(),
                norms: norms.clone(),
                faces: faces.clone(),
            }
        })
        .collect()
}

#[derive(Clone, Debug)]
pub struct Mesh {
    /// List of indexes pointing to the polygons in the list of shapes
    pub name: String,
    pub idx: PrimitiveIdx,
    pub polygons: Vec<usize>,
    pub verts: Rc<RefCell<Vec<Point3<f64>>>>,
    pub norms: Rc<RefCell<Vec<Dir3>>>,
    pub faces: Rc<RefCell<Vec<IdxTriangle>>>,
}

pub struct FaceIntersection {
    pub tri_u: PrimitiveIdx,
    pub tri_v: PrimitiveIdx,
    pub split_edges: SplitEdges,
}

pub struct MeshIntersection {
    pub mesh_u: PrimitiveIdx,
    pub mesh_v: PrimitiveIdx,
    pub tri_intersections: Vec<FaceIntersection>,
}

pub struct VertexMatch {
    pub meshes: (PrimitiveIdx, PrimitiveIdx),
    pub tris: (PrimitiveIdx, PrimitiveIdx),
    pub from: (IdxEdge, IdxEdge),
}

impl Collide<Mesh> for Mesh {
    fn overlap(&self, other: &Mesh) -> bool {
        return true;
        todo!()
    }
}

impl Split<Mesh, Vec<FaceIntersection>> for Mesh {
    type Inst = Mesh;
    fn intersect(&self, other: &Mesh) -> Vec<FaceIntersection> {
        // If there is no intersection, no point to proceed with checking pairwise triangles
        if !self.overlap(other) {
            return Vec::new();
        }

        //  Extract the triangles from mesh, keeping track of their idx for regrouping
        //  intersections
        let tri_u_union: Vec<_> = self
            .polygons
            .iter()
            .map(|idx| {
                (
                    PrimitiveIdx::Global(*idx),
                    self.faces.borrow()[*idx].clone(),
                )
            })
            .collect();
        let tri_v_union: Vec<_> = other
            .polygons
            .iter()
            .map(|idx| {
                (
                    PrimitiveIdx::Global(*idx),
                    other.faces.borrow()[*idx].clone(),
                )
            })
            .collect();

        // Find all intersections
        let mut intersections: Vec<FaceIntersection> = Vec::new();
        for (idx_tri_u, tri_u) in tri_u_union.iter() {
            for (idx_tri_v, tri_v) in tri_v_union.iter() {
                // FIXME: Would be better to have a collision tree checking here instead of
                // brute force pairwise collision tests
                // FIXME: Don't need to convert actually
                if tri_u.tri().overlap(&tri_v.tri()) {
                    let new_intersection = tri_u.intersect(&tri_v);
                    assert!(
                        !new_intersection.0.is_empty(),
                        "Intersection of overlapping triangles should not be empty"
                    );
                    intersections.push(FaceIntersection {
                        tri_u: *idx_tri_u,
                        tri_v: *idx_tri_v,
                        split_edges: new_intersection,
                    });
                }
            }
        }
        // TODO: Allocated new vertices with PrimitiveIdx::Local() to the global verts list
        // and rename the verts in order to avoid duplication in the case of same intersections
        // from distinct triangles (i.e. same edges but different triangles pair)
        intersections
    }

    fn split(&self, other: &Vec<FaceIntersection>) -> (Mesh, Vec<FaceIntersection>) {
        let tri_u_union: Vec<_> = self
            .polygons
            .iter()
            .map(|idx| {
                (
                    PrimitiveIdx::Global(*idx),
                    self.faces.borrow()[*idx].clone(),
                )
            })
            .collect();
        let mut new_mesh = Mesh {
            name: self.name.clone(),
            idx: self.idx,
            polygons: Vec::new(),
            verts: self.verts.clone(),
            norms: self.norms.clone(),
            faces: self.faces.clone(),
        };
        // NOTE: By this point all vertices in TriangleIntersection must be resolved and allocated
        // globally
        for (idx_tri_u, tri_u) in tri_u_union.iter() {
            let mut splits = SplitEdges(vec![]);
            other.iter().for_each(|intersection| {
                if intersection.tri_u == *idx_tri_u || intersection.tri_v == *idx_tri_u {
                    splits = splits.union(&intersection.split_edges);
                }
            });
            let new_tris_u = tri_u.split(&splits);

            //new_mesh.push(new_tris_u);
        }
        todo!()
    }
}
