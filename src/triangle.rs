use log::{trace, warn};
use nalgebra::{Point3, Vector3, Unit};

use crate::{Aabb, Collide,Segment};
type Dir3 = Unit<Vector3<f64>>;

#[derive(Debug, Clone)]
pub struct Triangle {
    pub verts: [Point3<f64>; 3],
    pub plane_norm: Dir3,
}

impl Triangle {

    pub fn new(verts: [Point3<f64>; 3]) -> Self {
        let plane_norm = Self::init_plane_norm(&verts);
        Self::new_with_norm(verts, plane_norm)
    }

    pub fn new_with_norm(verts: [Point3<f64>; 3], plane_norm: Dir3) -> Self {
        {
            let a = verts[0] - verts[1];
            let b = verts[1] - verts[2];
            assert!(
                a.cross(&b).abs().sum() > f64::EPSILON,
                "Edges of triangle can't be colinear: {:?}",
                verts
            );
        }
        Self { verts, plane_norm }
    }

    /// Initialise the plane normal.
    fn init_plane_norm(verts: &[Point3<f64>; 3]) -> Dir3 {
        Dir3::new_normalize((verts[0] - verts[2]).cross(&(verts[1] - verts[0])))
    }

    pub fn aabb(&self) -> Aabb {
        let mins = self.verts.iter().fold(
            Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            |acc, v| Point3::new(acc[0].min(v[0]), acc[1].min(v[1]), acc[2].min(v[2])),
        );
        let maxs = self.verts.iter().fold(
            Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            |acc, v| Point3::new(acc[0].max(v[0]), acc[1].max(v[1]), acc[2].max(v[2])),
        );
        Aabb{mins, maxs}
    }

    pub fn edges(&self) -> [Segment; 3] {
        [
            Segment::new(self.verts[1], self.verts[2]),
            Segment::new(self.verts[2], self.verts[0]),
            Segment::new(self.verts[0], self.verts[1]),
        ]
    }

    pub fn vertex_in(&self, vertex: Point3<f64>) -> bool {
        let segs = self.edges();
        let mut sign: Option<bool> = None;

        // 1. First check that vertex is in the plane of the triangle
        let vertex_to_tri = vertex - self.verts[0];
        if self.plane_norm.dot(&vertex_to_tri).abs() > 1e-9 {
            return false;
        }

        // 2. Then check that vertex is on the same side of all edges of the triangle
        for seg in &segs {
            let edge_vec = seg.end - seg.start;
            let to_vertex_vec = vertex - seg.start;
            let cross_prod = edge_vec.cross(&to_vertex_vec);
            // 2. Check if vertex is on an edge of the triangle
            if cross_prod.abs().sum() < 1e-9 {
                trace!("Vertex {:?} is on edge {:?} of triangle.", vertex, seg);
                return false;
            }
            // 3. Check that the vertex is on the same side of each segment
            let dir = self.plane_norm.dot(&cross_prod) > 0.0;
            match sign {
                None => {
                    sign = Some(dir);
                }
                Some(sign) => {
                    if sign != dir {
                        return false;
                    }
                }
            };
        }
        true
    }

}

impl Collide<Point3<f64>> for Triangle {
    // NOTE: Similar to Triangle.vertex_in, but including coincidence with triangle vertices and
    // edges
    fn overlap(&self, point: &Point3<f64>) -> bool {
        let segs = self.edges();
        let mut sign: Option<bool> = None;

        // 1. First check that vertex is in the plane of the triangle
        let vertex_to_tri = point - self.verts[0];
        if self.plane_norm.dot(&vertex_to_tri).abs() > 1e-9 {
            return false;
        }

        // 2. Then check that vertex is on the same side of all edges of the triangle
        for seg in &segs {
            let edge_vec = seg.end - seg.start;
            let to_vertex_vec = point - seg.start;
            if to_vertex_vec.abs().sum() < 1e-9 {
                // Vertex coincides with triangle vertex
                break;
            }
            let cross_prod = edge_vec.cross(&to_vertex_vec);
            // 3. Check that the vertex is on the same side of each segment
            let cross_prod_mag = self.plane_norm.dot(&cross_prod);
            if cross_prod_mag.abs() < 1e-9 {
                // Point is on the edge
                continue;
            }
            let dir = self.plane_norm.dot(&cross_prod) > 0.0;
            match sign {
                None => {
                    sign = Some(dir);
                }
                Some(sign) => {
                    if sign != dir {
                        return false;
                    }
                }
            };
        }
        true
    }
}

impl Collide<Triangle> for Triangle {
    fn overlap(&self, other: &Triangle) -> bool {
        const EPS_INTERSECTION: f64 = 1e-9;

        // 1. Check that triangle planes are parallel
        let planes_allignment = self.plane_norm.dot(&other.plane_norm);
        if planes_allignment.abs() < 0.999 {
            return false;
        }

        // 2. Check that the Aabb of the triangles has an intersection
        if !self.aabb().overlap(&other.aabb()) {
            return false;
        }

        // Choose maximal cross triangle distance
        let cross_triangle_edge = self
            .verts
            .iter()
            .zip(other.verts.iter())
            .map(|(u, v)| u - v)
            .max_by(|a, b| a.dot(a).partial_cmp(&b.dot(b)).unwrap())
            .unwrap();
        //let cross_triangle_edge = self.verts[ALPHA] - other.verts[ALPHA];

        // 3. Check that the triangles are coplanar
        if self.plane_norm.dot(&cross_triangle_edge).abs() > EPS_INTERSECTION {
            false
        } else {
            if planes_allignment > 0.0 {
                warn!(
                    "Triangles normals face the same way for: {:?} and {:?}",
                    self, other
                );
            }

            // 4. Check that either at least one vertex is inside the triangle or a cross edge
            //    intersection
            for &v in &other.verts {
                if self.vertex_in(v) {
                    trace!("Found vertex {:?} inside triangle", v);
                    return true;
                }
            }
            let segs_u = [
                Segment::new(other.verts[0], other.verts[1]),
                Segment::new(other.verts[1], other.verts[2]),
                Segment::new(other.verts[2], other.verts[0]),
            ];
            for seg_u in segs_u {
                // FIXME: This should exclude touching but not overlaping triangles
                if self.overlap(&seg_u) {
                    trace!("Found edge {:?} intersection with triangle", seg_u);
                    return true;
                }
            }
            false
        }
    }
}

impl Collide<Segment> for Triangle {
    fn overlap(&self, seg: &Segment) -> bool {
        let segs_u = self.edges();
        for seg_u in segs_u {
            // FIXME: Should be open intersection, however this seems to not detect any
            // intersection
            if let Some(_intersection) = seg.intersect(&seg_u) {
                return true;
            }
        }
        false
    }
}

impl Collide<Aabb> for Triangle {
    fn overlap(&self, cube: &Aabb) -> bool {
        let c = cube.centre();
        let e = cube.half_widths();

        let v0 = self.verts[0] - c;
        let v1 = self.verts[1] - c;
        let v2 = self.verts[2] - c;

        let f0 = v1 - v0;
        let f1 = v2 - v1;
        let f2 = v0 - v2;

        let u0 = Vector3::x_axis();
        let u1 = Vector3::y_axis();
        let u2 = Vector3::z_axis();

        let axis_test = |axis: &Vector3<f64>| {
            let p0 = v0.dot(axis);
            let p1 = v1.dot(axis);
            let p2 = v2.dot(axis);

            let r = e[2].mul_add(
                u2.dot(axis).abs(),
                e[0]
                    .mul_add(u0.dot(axis).abs(), e[1] * u1.dot(axis).abs()),
            );

            if (-(p0.max(p1).max(p2))).max(p0.min(p1).min(p2)) > r {
                return false;
            }

            true
        };

        if !axis_test(&u0) {
            return false;
        }
        if !axis_test(&u1) {
            return false;
        }
        if !axis_test(&u2) {
            return false;
        }

        let axis_u0_f0 = u0.cross(&f0);
        let axis_u0_f1 = u0.cross(&f1);
        let axis_u0_f2 = u0.cross(&f2);

        let axis_u1_f0 = u1.cross(&f0);
        let axis_u1_f1 = u1.cross(&f1);
        let axis_u1_f2 = u1.cross(&f2);

        let axis_u2_f0 = u2.cross(&f0);
        let axis_u2_f1 = u2.cross(&f1);
        let axis_u2_f2 = u2.cross(&f2);

        if !axis_test(&axis_u0_f0) {
            return false;
        }
        if !axis_test(&axis_u0_f1) {
            return false;
        }
        if !axis_test(&axis_u0_f2) {
            return false;
        }

        if !axis_test(&axis_u1_f0) {
            return false;
        }
        if !axis_test(&axis_u1_f1) {
            return false;
        }
        if !axis_test(&axis_u1_f2) {
            return false;
        }

        if !axis_test(&axis_u2_f0) {
            return false;
        }
        if !axis_test(&axis_u2_f1) {
            return false;
        }
        if !axis_test(&axis_u2_f2) {
            return false;
        }

        if !axis_test(&self.plane_norm) {
            return false;
        }

        true
    }
}
