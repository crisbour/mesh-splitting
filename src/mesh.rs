use std::{cell::RefCell, collections::{BTreeMap, HashMap}, ops::Deref, rc::Rc};

use log::info;
use nalgebra::{Point3, Unit, Vector3};
use obj::Obj;

use crate::{
    Collide, IdxTriangle, Split, SplitEdges,
    primitives::{IdxEdge, IdxIntersection, Normal, PrimitiveIdx, Vertex},
};

type Dir3 = Unit<Vector3<f64>>;

#[derive(Clone, Debug)]
pub struct Verts {
    verts: Rc<RefCell<Vec<Point3<f64>>>>,
}

impl Deref for Verts {
    type Target = Rc<RefCell<Vec<Point3<f64>>>>;
    fn deref(&self) -> &Rc<RefCell<Vec<Point3<f64>>>> {
        &self.verts
    }
}

impl Verts {
    pub fn new() -> Self {
        let verts = Rc::new(RefCell::new(Vec::new()));
        Self { verts }
    }
    pub fn with_capacity(capacity: usize) -> Self {
        let verts = Rc::new(RefCell::new(Vec::with_capacity(capacity)));
        Self { verts }
    }
    pub fn alllocate_size(&self, size: usize) {
        let current_size = self.verts.borrow().len();
        self.verts.borrow_mut().resize(current_size + size, Point3::new(0.0, 0.0, 0.0));
    }
    pub fn allocate(&self, vert: Point3<f64>) -> PrimitiveIdx {
        let idx = self.verts.borrow().len();
        self.verts.borrow_mut().push(vert);
        PrimitiveIdx::Global(idx)
    }
    pub fn allocate_iter(&self, verts: impl IntoIterator<Item=Point3<f64>>) -> Vec<PrimitiveIdx> {
        let base_idx = self.verts.borrow().len();
        verts.into_iter().enumerate().map(|(i, vert)| {
            self.verts.borrow_mut().push(vert);
            PrimitiveIdx::Global(base_idx + i)
        }).collect()
    }
}

#[derive(Clone, Debug)]
pub struct Norms {
    norms: Rc<RefCell<Vec<Dir3>>>,
}

impl Deref for Norms {
    type Target = Rc<RefCell<Vec<Dir3>>>;
    fn deref(&self) -> &Rc<RefCell<Vec<Dir3>>> {
        &self.norms
    }
}

impl Norms {
    pub fn new() -> Self {
        let norms = Rc::new(RefCell::new(Vec::new()));
        Self { norms }
    }
    pub fn with_capacity(capacity: usize) -> Self {
        let norms = Rc::new(RefCell::new(Vec::with_capacity(capacity)));
        Self { norms }
    }
    pub fn alllocate_size(&self, size: usize) {
        let current_size = self.norms.borrow().len();
        self.norms.borrow_mut().resize(current_size + size, Unit::new_normalize(Vector3::new(1.0, 0.0, 0.0)));
    }
    pub fn allocate(&self, dir: Dir3) -> PrimitiveIdx {
        let idx = self.norms.borrow().len();
        self.norms.borrow_mut().push(dir);
        PrimitiveIdx::Global(idx)
    }
    pub fn allocate_iter(&self, dirs: impl IntoIterator<Item=Dir3>) -> Vec<PrimitiveIdx> {
        let base_idx = self.norms.borrow().len();
        dirs.into_iter().enumerate().map(|(i, dir)| {
            self.norms.borrow_mut().push(dir);
            PrimitiveIdx::Global(base_idx + i)
        }).collect()
    }
}

#[derive(Clone, Debug)]
pub struct Faces {
    faces: Rc<RefCell<Vec<IdxTriangle>>>,
}

impl Deref for Faces {
    type Target = Rc<RefCell<Vec<IdxTriangle>>>;
    fn deref(&self) -> &Rc<RefCell<Vec<IdxTriangle>>> {
        &self.faces
    }
}

impl Faces {
    pub fn new() -> Self {
        let faces = Rc::new(RefCell::new(Vec::new()));
        Self { faces }
    }
    pub fn with_capacity(capacity: usize) -> Self {
        let faces = Rc::new(RefCell::new(Vec::with_capacity(capacity)));
        Self { faces }
    }
    // FIXME: Allocate size should increase capacity, not resize and initialize to default
    pub fn alllocate_size(&self, size: usize) {
        let current_size = self.faces.borrow().len();
        self.faces.borrow_mut().resize(current_size + size, IdxTriangle::default());
    }
    pub fn allocate(&self, tri: IdxTriangle) -> PrimitiveIdx {
        let idx = self.faces.borrow().len();
        self.faces.borrow_mut().push(tri);
        PrimitiveIdx::Global(idx)
    }
    pub fn allocate_iter(&self, tris: impl IntoIterator<Item=IdxTriangle>) -> Vec<PrimitiveIdx> {
        let base_idx = self.faces.borrow().len();
        tris.into_iter().enumerate().map(|(i, tri)| {
            self.faces.borrow_mut().push(tri);
            PrimitiveIdx::Global(base_idx + i)
        }).collect()
    }
}

pub fn parse_obj(
    obj: Obj,
) -> (Vec<Mesh>, Verts, Norms, Faces) {
    let verts = Verts::with_capacity(obj.data.position.len());
    let norms = Norms::with_capacity(obj.data.normal.len());
    let faces = Faces::new();

    let allocated_verts_idx = verts.allocate_iter(
        obj
        .data
        .position
        .iter()
        .map(|vs| {
            let vs_f64: [f64; 3] = vs.map(|v| v as f64);
            Point3::new(vs_f64[0], vs_f64[1], vs_f64[2])
        })
    );
    println!("Allocated verts idx: {:?}", allocated_verts_idx);

    let allocated_norms_idx = norms.allocate_iter(
        obj
        .data
        .normal
        .iter()
        .map(|vs| {
            let vs_f64: [f64; 3] = vs.map(|v| v as f64);
            Dir3::new_normalize(Vector3::new(vs_f64[0], vs_f64[1], vs_f64[2]))
        })
    );

    println!("Allocated norms idx: {:?}", allocated_norms_idx);

    let meshes =
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

            //faces.alllocate_size(obj_faces.len());
            let polygons = faces.allocate_iter(obj_faces);
            println!("Allocated faces idx for mesh {}: {:?}", obj.name, polygons);

            Mesh {
                name: obj.name.clone(),
                idx: PrimitiveIdx::Global(idx),
                polygons,
                verts: verts.clone(),
                norms: norms.clone(),
                faces: faces.clone(),
            }
        })
        .collect();
    (meshes, verts, norms, faces)
}

#[derive(Clone, Debug)]
pub struct Mesh {
    /// List of indexes pointing to the polygons in the list of shapes
    pub name: String,
    pub idx: PrimitiveIdx,
    pub polygons: Vec<PrimitiveIdx>,
    pub verts: Verts,
    pub norms: Norms,
    pub faces: Faces,
}

impl Mesh {
    pub fn push(&mut self, new_tris: Vec<IdxTriangle>) {
        let new_polygons = self.faces.allocate_iter(new_tris);
        self.polygons.extend(new_polygons);
    }
}

#[derive(Clone)]
pub struct FaceIntersection {
    pub tri_u: PrimitiveIdx,
    pub tri_v: PrimitiveIdx,
    pub split_edges: SplitEdges,
}

impl FaceIntersection {
    pub fn rename_vertices(&mut self, idx_alloc_map: &HashMap<IdxIntersection, PrimitiveIdx>) {
        self.split_edges.rename_vertices(idx_alloc_map);
    }
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
            .map(|&idx| {
                (
                    idx,
                    self.faces.borrow()[*idx].clone(),
                )
            })
            .collect();
        let tri_v_union: Vec<_> = other
            .polygons
            .iter()
            .map(|&idx| {
                (
                    idx,
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
                    info!("Intersect triangles {:?} and {:?} from meshes {:?} and {:?}",
                        idx_tri_u, idx_tri_v, self.idx, other.idx);
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
        intersections
    }

    fn split(&self, other: Vec<FaceIntersection>) -> (Mesh, Vec<FaceIntersection>) {
        let tri_u_union: Vec<_> = self
            .polygons
            .iter()
            .map(|&idx| {
                (
                    idx,
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

        // Map from local vertex idx in intersection to new global vertex idx in new_mesh
        let mut idx_alloc_map: HashMap<IdxIntersection, PrimitiveIdx> = HashMap::new();

        // Collect split edges derived from the target triangle
        for (idx_tri_u, tri_u) in tri_u_union.iter() {
            // Collect all splits with tri_u, identified by idx_tri_u:PrimitiveIdx
            let mut splits = SplitEdges(vec![]);
            // TODO: Instead of iterating through all, sort intersection and binary search for the
            // ones matching idx_tri_u
            other.iter().for_each(|intersection| {
                if intersection.tri_u == *idx_tri_u || intersection.tri_v == *idx_tri_u {
                    splits = splits.union(&intersection.split_edges);
                }
            });

            // Allocate vertices from splits, and keep track of the mapping from local vertex idx to global
            let local_verts = Vertex::from_edges(&splits.0)
                .into_iter()
                .filter(|vert| matches!(vert.idx, PrimitiveIdx::Local(..)));

            for vert in local_verts {
                assert!(vert.from.is_some(), "Local vertex must have intersection edges described");
                if !idx_alloc_map.contains_key(&vert.from.unwrap()) {
                    let new_idx = new_mesh.verts.allocate(vert.value);
                    debug_assert!(matches!(new_idx, PrimitiveIdx::Global(_)), "Vertices from edges should be allocated globally");
                    idx_alloc_map.insert(vert.from.unwrap(), new_idx);
                }
            }

            // Rneame vertices idx, mapping to the allocated global idx
            splits.rename_vertices(&idx_alloc_map);

            // Sanity check: All vertices in splits should now be allocated globally
            Vertex::from_edges(&splits.0)
                .into_iter()
                .for_each(|vert| {
                    assert!(matches!(vert.idx, PrimitiveIdx::Global(_)), "All vertices in splits should be allocated globally");
                });

            // Split the triangle with the splits, and push to the new mesh
            let (new_tris_u, _splits) = tri_u.split(splits);
            new_mesh.push(new_tris_u);
        }

        let new_face_inters = other
            .iter()
            .map(|face_inter| {
                let mut new_face_inter = face_inter.clone();
                new_face_inter.rename_vertices(&idx_alloc_map);
                new_face_inter
            })
            .collect();

        (new_mesh, new_face_inters)
    }
}
