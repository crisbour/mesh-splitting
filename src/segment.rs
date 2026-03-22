use core::f64;

use crate::{
    Collide, IdxTriangle, Triangle,
    primitives::{Edge, IdxIntersection, Polygon, PrimitiveIdx, Vertex},
};
use anyhow::{Context, Result};
use log::{error, trace};
use nalgebra::Point3;

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
        let e = self.end - self.start;
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

        let det = a * c - b * b;
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
            //trace!("Segments do not intersect within their lengths. alpha_self: {}, alpha_other: {}", alpha_u, alpha_v);
            return None;
        }

        let p_u = self.at(alpha_u);
        let p_v = other.at(alpha_v);

        let p_diff = p_u - p_v;
        // Check that closest points actually coincide
        if p_diff.dot(&p_diff) <= EPS_INTERSECT * EPS_INTERSECT {
            Some(p_u)
        } else {
            trace!(
                "Closest points do not coincide. p_self: {:?}, p_other: {:?}, p_diff: {:?}",
                p_u, p_v, p_diff
            );
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
            //trace!("Segments do not intersect within their lengths. alpha_self: {}, alpha_other: {}", alpha_u, alpha_v);
            return None;
        }

        let fuzzy = alpha_v < eps || alpha_v > 1.0 - eps;

        let p_u = self.at(alpha_u);
        let p_v = other.at(alpha_v);

        let p_diff = p_u - p_v;
        // Check that closest points actually coincide
        if p_diff.dot(&p_diff) <= eps * eps {
            Some((p_u, fuzzy))
        } else {
            trace!(
                "Closest points do not coincide. p_self: {:?}, p_other: {:?}, p_diff: {:?}",
                p_u, p_v, p_diff
            );
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

impl From<Edge> for Segment {
    fn from(edge: Edge) -> Self {
        Segment::new(edge.0.value, edge.1.value)
    }
}

impl Collide<Point3<f64>> for Segment {
    // FIXME: Perhaps increase EPS to allow for some numerical precision issues
    #[inline]
    fn overlap(&self, other: &Point3<f64>) -> bool {
        let e = self.end - self.start;
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

#[derive(Clone, Debug)]
pub struct SplitEdges {
    /// All edges produced by intersection, including those after triangulation of the polygon(s)
    pub edges: Vec<Edge>,
    /// Perimeter edges of the intersection(s) , as indices into `edges`
    pub outer: Vec<usize>,
}

impl SplitEdges {
    /// Construct a polygon from edges
    /// TODO: Check that edges form a valid polygon, and that the order of vertices is consistent (e.g. counter-clockwise)
    pub fn new(edges: Vec<Edge>) -> Self {
        let outer = (0..edges.len()).collect();
        Self { edges, outer }
    }
    /// Return perimeter edges
    pub fn outer_edges(&self) -> Vec<Edge> {
        self.outer.iter().map(|&i| self.edges[i].clone()).collect()
    }
    /// Triangulate the polygon formed by the outer edges, and add the resulting segments to the SplitEdges
    ///   - WARN: This works only if outer edges describe a convex polygon
    pub fn triangulate(&self) -> Result<Self> {
        let convex_polygon = Polygon::try_from(self.outer_edges().clone())
            .context("Outer edges must form a valid polygon for triangulation")?;

        let triangles = convex_polygon.triangulate();

        let mut edges = self.edges.clone();
        for tri in triangles {
            let tri_edges = [
                Edge::new(tri.verts[0], tri.verts[1]),
                Edge::new(tri.verts[1], tri.verts[2]),
                Edge::new(tri.verts[2], tri.verts[0]),
            ];
            for edge in tri_edges {
                if !edges.contains(&edge) {
                    edges.push(edge);
                }
            }
        }
        Ok(Self {
            edges,
            outer: self.outer.clone(),
        })
    }
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }
    /// Collect all edges from both SplitEdges sets,
    /// renaming (allocate PrimitiveIdx::Local freshly) vertices produced from intersection
    /// such that we match segments.
    ///   - Outer edges are updated to include only those that are present in one of the sets, but
    ///   not both, hence producing boundaries even for unions of disjoint sets
    ///   - WARN: The sets to be united must be disjoint in the areas they cover
    pub fn union(&self, other: &SplitEdges) -> SplitEdges {
        trace!("Union of SplitEdges: {:?} and {:?}", self, other);
        let mut edges = self.edges.clone();
        let mut verts = Vertex::from_edges(&edges);
        let mut outer = self.outer.clone();

        // 1. Find starting local index
        let mut local_idx = 0;
        for edge in self.edges.iter() {
            match edge.0.idx {
                PrimitiveIdx::Local(idx) => {
                    if local_idx <= idx {
                        local_idx = idx+1;
                    }
                }
                _ => {}
            }
            match edge.1.idx {
                PrimitiveIdx::Local(idx) => {
                    if local_idx <= idx {
                        local_idx = idx+1;
                    }
                }
                _ => {}
            }
        }

        // 2. Vertex renaming for set to be merged in current one
        let mut other_edges = other.edges.clone();
        for other_edge in other_edges.iter_mut() {
            // Rename PrimitiveIdx::Local, but preserve Global
            for other_vert in [&mut other_edge.0, &mut other_edge.1] {
                if let PrimitiveIdx::Local(_idx) = other_vert.idx {
                    if let Some(vert) = verts.iter().find(|&v| v == other_vert) {
                        if let PrimitiveIdx::Global(_) = other_vert.idx {
                            assert_eq!(other_vert.idx, vert.idx, "Expected {:?} not to be renamed", other_vert.idx);
                        }
                        trace!(
                            "{:?} already exists,\
                             reusing existing vertex {:?} instead of creating a new one.",
                            other_vert,
                            vert
                        );
                        *other_vert = *vert;
                    } else {
                        let new_vert = Vertex {
                            value: other_vert.value,
                            idx: PrimitiveIdx::Local(local_idx),
                            from: other_vert.from,
                        };
                        trace!(
                            "Renaming vertex {:?} to {:?} for union, as it is not present in the current set",
                            other_vert.idx, local_idx
                        );
                        local_idx += 1;
                        verts.push(new_vert);
                        *other_vert = new_vert;
                    }
                }
            }
        }

        // 3. Merge edges, and update outer edges
        for (other_pos, other_edge) in other_edges.into_iter().enumerate() {
            let pos = edges.iter().position(|&e| e == other_edge);
            if let Some(pos) = pos {
                if let Some(outer_pos) = outer.iter().position(|&p| p == pos) {
                    trace!(
                        "Removing edge {:?} from outer edges, as it is shared between both SplitEdges sets.",
                        other_edge
                    );
                    outer.remove(outer_pos);
                } else if other.outer.contains(&other_pos) {
                    error!("Did not expect outer edge of the other set to be inside the merging one.\n
                        That means the set are not disjoint in surface");
                }
            } else {
                if other.outer.contains(&other_pos){
                    outer.push(edges.len());
                }
                edges.push(other_edge);
            }
        }

        // Make sure we have a valid union
        assert!(
            outer.len() <= edges.len(),
            "Outer edges cannot be more than total edges"
        );
        assert!(
            outer.len() >= 3,
            "Outer edges must be at least 3 to form a polygon, but found: {:?}",
            outer
        );
        SplitEdges { edges, outer }
    }

    pub fn rename_vertices(
        &mut self,
        vertex_map: &std::collections::HashMap<IdxIntersection, PrimitiveIdx>,
    ) {
        for edge in &mut self.edges {
            match edge.0.from {
                Some(idx_inter) => {
                    if let Some(&new_idx) = vertex_map.get(&idx_inter) {
                        edge.0.idx = new_idx;
                    }
                }
                _ => {}
            }
            match edge.1.from {
                Some(idx_inter) => {
                    if let Some(&new_idx) = vertex_map.get(&idx_inter) {
                        edge.1.idx = new_idx;
                    }
                }
                _ => {}
            }
        }
    }
}

impl From<Vec<IdxTriangle>> for SplitEdges {
    fn from(value: Vec<IdxTriangle>) -> Self {
        let mut edges = Vec::new();
        for tri in value {
            let tri_edges = [
                Edge::new(tri.verts[0], tri.verts[1]),
                Edge::new(tri.verts[1], tri.verts[2]),
                Edge::new(tri.verts[2], tri.verts[0]),
            ];
            for edge in tri_edges {
                if !edges.contains(&edge) {
                    edges.push(edge);
                }
            }
        }
        SplitEdges::new(edges)
    }
}

impl Collide<SplitEdges> for IdxTriangle {
    fn overlap(&self, other: &SplitEdges) -> bool {
        let split_verts = Vertex::from_edges(&other.edges);
        let tri: Triangle = self.clone().into();
        for vert in split_verts {
            if tri.vertex_in(vert.value) {
                return true;
            }
        }
        other.edges.iter().any(|edge| self.overlap(edge))
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

    #[test]
    fn test_segment_intersect() {
        let seg1 = Segment::new(Point3::new(0., 0., 0.), Point3::new(1., 1., 0.));
        let seg2 = Segment::new(Point3::new(0., 1., 0.), Point3::new(1., 0., 0.));
        let intersection = seg1.intersect(&seg2);
        assert!(intersection.is_some());
        assert!((intersection.unwrap()-Point3::new(0.5, 0.5, 0.)).abs().sum() < 1e-9,
            "Expected intersection at (0.5, 0.5, 0), got {:?}",
            intersection.unwrap());

        let seg1 = Segment::new(Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0));
        let seg2 = Segment::new(Point3::new(0.8, 0.5, 0.0), Point3::new(2.0, 0.0, 0.0));
        assert!(seg1.intersect(&seg2).is_some());
    }

    #[test]
    fn test_segment_intersect_open() {
        let seg1 = Segment::new(Point3::new(0., 0., 0.), Point3::new(1., 1., 0.));
        let seg2 = Segment::new(Point3::new(0., 1., 0.), Point3::new(1., 0., 0.));
        let intersection = seg1.intersect_open(&seg2);
        assert!(intersection.is_some());
        assert!((intersection.unwrap()-Point3::new(0.5, 0.5, 0.)).abs().sum() < 1e-9,
            "Expected intersection at (0.5, 0.5, 0), got {:?}",
            intersection.unwrap());
    }

    #[test]
    fn test_segment_non_intersect_open() {
        let seg1 = Segment::new(Point3::new(0., 0., 0.), Point3::new(1., 1., 0.));
        let seg2 = Segment::new(Point3::new(1., 1., 0.), Point3::new(2., 0., 0.));
        let intersection = seg1.intersect_open(&seg2);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_segment_non_intersect() {
        let seg1 = Segment::new(Point3::new(0., 0., 0.), Point3::new(1., 1., 0.));
        let seg2 = Segment::new(Point3::new(0., 1., 0.), Point3::new(1., 2., 0.));
        let intersection = seg1.intersect(&seg2);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_union_split_segments () {
        // TODO: Write tests of interest, but also edge cases for degenerate scenarious

    }
}
