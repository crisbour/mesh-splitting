use core::f64;

use log::trace;
use nalgebra::Point3;
use crate::{Collide, IdxTriangle, Triangle, primitives::{Edge, PrimitiveIdx, Vertex}};

#[derive(Debug)]
pub struct Segment {
    pub start: Point3<f64>,
    pub end: Point3<f64>,
}

impl Segment {
    pub fn new(start: Point3<f64>, end: Point3<f64>) -> Self {
        Self { start, end }
    }

    pub fn length(&self) -> f64 {
        nalgebra::distance(&self.start, &self.end)
    }

    #[inline]
    pub fn at(&self, alpha: f64) -> Point3<f64> {
        self.start + (self.end - self.start) * alpha
    }

    pub fn parametric_dist(&self, p: &Point3<f64>) -> f64 {
        let e = self.end -self.start;
        (p - self.start).dot(&e) / e.dot(&e) // alpha
    }

    #[inline]
    fn intersect_unchecked(&self, other: &Segment) -> Option<(f64, f64)> {
        const EPS_INTERSECT: f64 = 1e-9;

        // TODO:Check this algorithms matches with the same method used for Barycentric solution

        // Check if segments are colinear
        let u = self.end - self.start;
        let v = other.end - other.start;

        if u.cross(&v).abs().sum() < f64::EPSILON { // norm1
            return None;
        }

        let u = self.end - self.start;
        let v = other.end - other.start;
        let w0 = self.start - other.start;

        let a = u.dot(&u); // u · u
        let b = u.dot(&v); // u · v
        let c = v.dot(&v); // v · v
        let d = u.dot(&w0); // u · (u1-v1)
        let e = v.dot(&w0); // v · (u1-v1)

        let det =  a * c - b * b;
        if det.abs() < EPS_INTERSECT {
            return None;
        }

        let alpha_u = (b * e - c * d) / det;
        let alpha_v = (a * e - b * d) / det;
        Some((alpha_u, alpha_v))
    }

    /// Calculate the intersection point of two line segments in 3D space.
    pub fn intersect(&self, other: &Segment) -> Option<Point3<f64>> {
        const EPS_INTERSECT: f64 = 1e-9;

        // Check if segments are colinear
        let u = self.end - self.start;
        let v = other.end - other.start;
        if u.cross(&v).abs().sum() < f64::EPSILON { // norm1
            return None;
        }

        let (alpha_u, alpha_v) = self.intersect_unchecked(&other)?;

        if alpha_u < 0.0 || alpha_u > 1.0 || alpha_v < 0.0 || alpha_v > 1.0 {
            trace!("Segments do not intersect within their lengths. alpha_self: {}, alpha_other: {}", alpha_u, alpha_v);
            return None;
        }

        let p_u = self.at(alpha_u);
        let p_v = other.at(alpha_v);

        let p_diff = p_u - p_v;
        // Check that closest points actually coincide
        if p_diff.dot(&p_diff) <= EPS_INTERSECT * EPS_INTERSECT {
            Some(p_u)
        } else {
            trace!("Closest points do not coincide. p_self: {:?}, p_other: {:?}, p_diff: {:?}", p_u, p_v, p_diff);
            None
        }
    }

    /// Calculate the intersection point of two line segments in 3D space.
    pub fn intersect_with_eps(&self, other: &Segment, eps: f64) -> Option<(Point3<f64>, bool)> {
        // Check if segments are colinear
        let u = self.end - self.start;
        let v = other.end - other.start;
        if u.cross(&v).abs().sum() < f64::EPSILON { // norm1
            return None;
        }

        let (alpha_u, alpha_v) = self.intersect_unchecked(&other)?;

        if alpha_u < 0.0 || alpha_u > 1.0 || alpha_v < -eps || alpha_v > 1.0 + eps {
            trace!("Segments do not intersect within their lengths. alpha_self: {}, alpha_other: {}", alpha_u, alpha_v);
            return None;
        }

        let fuzzy = alpha_v < eps || alpha_v > 1.0 - eps;

        let p_u = self.at(alpha_u);
        let p_v = other.at(alpha_v);

        let p_diff = p_u - p_v;
        // Check that closest points actually coincide
        if p_diff.dot(&p_diff) <= eps*eps {
            Some((p_u, fuzzy))
        } else {
            trace!("Closest points do not coincide. p_self: {:?}, p_other: {:?}, p_diff: {:?}", p_u, p_v, p_diff);
            None
        }
    }

    pub fn intersect_open(&self, other: &Segment) -> Option<Point3<f64>> {
        const EPS_INTERSECT: f64 = 1e-9;

        let (alpha_u, alpha_v) = self.intersect_unchecked(&other)?;

        if alpha_u < 0.0 || alpha_u > 1.0 || alpha_v <= 0.0 || alpha_v >= 1.0 {
            return None;
        }

        let p_u = self.at(alpha_u);
        let p_v = other.at(alpha_v);

        let p_diff = p_u - p_v;
        // Check that closest points actually coincide
        if p_diff.dot(&p_diff) <= EPS_INTERSECT * EPS_INTERSECT {
            Some(p_u)
        } else {
            None
        }
    }

    pub fn colinear(&self, other: &Segment) -> bool {
        let u = self.end - self.start;
        let v = other.end - other.start;
        let delta = other.end - self.start;
        u.cross(&v).abs().sum() <= f64::EPSILON && u.cross(&delta).abs().sum() <= f64::EPSILON
    }

    pub fn colinear_with_eps(&self, other: &Segment, eps: f64) -> bool {
        let u = self.end - self.start;
        let v = other.end - other.start;
        u.cross(&v).abs().sum() <= eps
    }
}

impl Collide<Point3<f64>> for Segment {
    // FIXME: Perhaps increase EPS to allow for some numerical precision issues
    #[inline]
    fn overlap(&self, other: &Point3<f64>) -> bool {
        let e = self.end -self.start;
        let alpha = (other - self.start).dot(&e) / e.dot(&e);
        let closest = self.at(alpha);
        let diff = closest - other;
        if diff.dot(&diff) < f64::EPSILON {
            alpha >= -f64::EPSILON && alpha <= 1.0 + f64::EPSILON
        } else {
            false
        }
    }
}

// ------------------------------------------------------------------------------------
// SplitEdges produced by polygon intersections after triangulation
// ------------------------------------------------------------------------------------

pub struct SplitEdges(pub Vec<Edge>);

impl SplitEdges {
    /// Extract the segments that are overlapping the target triangle
    pub fn intersect(&self, idx_tri: &IdxTriangle) -> SplitEdges {
        trace!("Splitting triangle with segments: {:?}", self.0);

        // Convert to simple spatial triangle
        let tri = idx_tri.tri();

        let mut verts: Vec<Vertex> = Vec::new();
        self.0.iter().for_each(|Edge(v0, v1)| {
            if !verts.contains(v0) && tri.overlap(&v0.value) {
                verts.push(*v0);
            }
            if !verts.contains(v1) && tri.overlap(&v0.value) {
                verts.push(*v1);
            }
        });

        // Only include segments that have both ends overlapping the triangle
        let edges = self.0
            .iter()
            .filter(|Edge(v0, v1)|
                tri.overlap(&v0.value) && tri.overlap(&v1.value)
            )
            .cloned()
            .collect();
        trace!("Segments after intersection with triangle: {:?}", edges);

        SplitEdges(edges)
    }
    /// Collect all edges from both SplitEdges sets,
    /// renaming (allocate PrimitiveIdx::Local freshly) vertices produced from intersection
    /// such that we match segments.
    pub fn union(&self, other: &SplitEdges) -> SplitEdges {
        let mut edges = self.0.clone();
        let mut local_idx = 0;
        for edge in self.0.iter() {
            match edge.0.idx {
                PrimitiveIdx::Local(idx) => {
                    if local_idx < idx {
                        local_idx = idx;
                    }
                },
                _ => {}
            }
            match edge.1.idx {
                PrimitiveIdx::Local(idx) => {
                    if local_idx < idx {
                        local_idx = idx;
                    }
                },
                _ => {}
            }
        }
        for other_edge in &other.0 {
            let pos = edges.iter().position(|e| e == other_edge);
            if let Some(pos) = pos {
                // Edge already exists, but check that the other vertices are not higher priority
                // i.e. self.edges vertex in Local, while the one from other.edges is Global
                let mut replace = false;
                let v0 = match other_edge.0.idx {
                    PrimitiveIdx::Global(_idx) => {
                        if matches!(edges[pos].0.idx, PrimitiveIdx::Local(_)) {
                            trace!("Replacing edge vertex {:?} with higher priority global vertex {:?}",
                                edges[pos].0, other_edge.0);
                            replace = true;
                        }
                        other_edge.0
                    },
                    _ => {
                        edges[pos].0
                    }
                };
                let v1 = match other_edge.1.idx {
                    PrimitiveIdx::Global(_idx) => {
                        if matches!(edges[pos].1.idx, PrimitiveIdx::Local(_)) {
                            trace!("Replacing edge vertex {:?} with higher priority global vertex {:?}",
                                edges[pos].1, other_edge.1);
                            replace = true;
                        }
                        other_edge.1
                    },
                    _ => {
                        edges[pos].1
                    }
                };
                if replace {
                    edges[pos] = Edge(v0, v1);
                }
            } else {
                // Rename PrimitiveIdx::Local, but preserve Global
                let v0 = match other_edge.0.idx {
                    PrimitiveIdx::Local(_idx) => {
                        local_idx += 1;
                        Vertex{
                            value: other_edge.0.value,
                            idx: PrimitiveIdx::Local(local_idx),
                            from: other_edge.0.from,
                        }
                    },
                    _ => {other_edge.0}
                };
                let v1 = match other_edge.1.idx {
                    PrimitiveIdx::Local(_idx) => {
                        local_idx += 1;
                        Vertex{
                            value: other_edge.1.value,
                            idx: PrimitiveIdx::Local(local_idx),
                            from: other_edge.1.from,
                        }
                    },
                    _ => {other_edge.1}
                };
                edges.push(Edge(v0, v1));
            }
            if !edges.contains(other_edge) {
            }
        }
        SplitEdges(edges)
    }
}

impl From<Vec<IdxTriangle>> for SplitEdges {
    fn from(value: Vec<IdxTriangle>) -> Self {
        let mut edges = Vec::new();
        for tri in value {
            let tri_edges = [Edge(tri.verts[0], tri.verts[1]), Edge(tri.verts[1], tri.verts[2]), Edge(tri.verts[2], tri.verts[0])];
            for edge in tri_edges {
                if !edges.contains(&edge) {
                    edges.push(edge);
                }
            }
        }
        SplitEdges(edges)
    }
}

impl From<&SplitEdges> for Vec<IdxTriangle> {
    fn from(value: &SplitEdges) -> Self {
        todo!()
    }
}

impl From<SplitEdges> for Vec<IdxTriangle> {
    fn from(value: SplitEdges) -> Self {
        todo!()
    }
}

impl Collide<SplitEdges> for IdxTriangle {
    fn overlap(&self, other: &SplitEdges) -> bool {
        let split_verts = Vertex::from_edges(&other.0);
        let tri: Triangle = self.clone().into();
        for vert in split_verts {
            if tri.vertex_in(vert.value) {
                return true;
            }
        }
        other
            .0
            .iter()
            .any(|edge| self.overlap(edge))
    }
}


#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    #[test]
    fn test_intersection() {
        let seg1 = Segment::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
        let seg2 = Segment::new(Point3::new(0.0, -1.0, 0.0), Point3::new(0.0, 1.0, 0.0));
        assert!(matches!(seg1.intersect(&seg2), Some(_)));
        assert!(matches!(seg2.intersect_open(&seg1), None));

        let seg2 = Segment::new(Point3::new(0.0, -1.0, 0.0), Point3::new(0.0, 0.0, 0.0));
        assert!(matches!(seg1.intersect(&seg2), Some(_)));
        assert!(matches!(seg1.intersect_open(&seg2), None));

        let seg2 = Segment::new(Point3::new(1e-9, -1.0, 0.0), Point3::new(1e-9, 1.0, 0.0));
        assert!(matches!(seg1.intersect(&seg2), Some(_)));
        assert!(matches!(seg2.intersect_open(&seg1), Some(_)));

    }

    #[test]
    fn test_colinear() {
        let seg1 = Segment::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
        let seg2 = Segment::new(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0));
        assert!(seg1.colinear(&seg2));
        assert!(seg1.colinear_with_eps(&seg2, 1e-9));
    }

    #[test]
    fn test_intersection_unchecked() {
        let seg1 = Segment::new(Point3::new(0.0, 1.0, 0.0), Point3::new(0.0, 0.0, 1.0));
        let seg2 = Segment::new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 2.0, 2.0));

        let (alpha_u, alpha_v) = seg1.intersect_unchecked(&seg2).unwrap();
        assert_approx_eq!(alpha_u, 0.5);
        assert_approx_eq!(alpha_v, 0.25);
        assert!(matches!(seg1.intersect(&seg2), Some(_)));
        assert!(matches!(seg1.intersect_open(&seg2), Some(_)));
        assert!(matches!(seg1.intersect_with_eps(&seg2, 1e-9), Some(_)));

    }

}

