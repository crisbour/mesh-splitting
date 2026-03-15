use std::{hash::Hash, ops::Deref};

use nalgebra::{Point3, Unit, Vector3};

use crate::IdxTriangle;

type Dir3 = Unit<Vector3<f64>>;

// Edge pointing to the index of the vertices in the list of vertices
#[derive(Debug, Clone, Copy)]
pub struct Edge(pub Vertex, pub Vertex);

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        (self.0 == other.0 && self.1 == other.1) || (self.0 == other.1 && self.1 == other.0)
    }
}
impl Eq for Edge {}

impl Edge {
    pub fn from_idx(idx_edge: IdxEdge, verts: &Vec<Point3<f64>>) -> Self {
        Edge(
            Vertex::new(*idx_edge.0, verts),
            Vertex::new(*idx_edge.1, verts),
        )
    }
}

#[derive(Clone, Debug, Copy, Hash)]
pub struct IdxEdge(pub PrimitiveIdx, pub PrimitiveIdx);

impl IdxEdge {
    pub fn new(idx1: PrimitiveIdx, idx2: PrimitiveIdx) -> Self {
        let min_idx = idx1.min(idx2);
        let max_idx = idx1.max(idx2);
        //if min_idx == max_idx {
        //    panic!("Edge cannot have the same vertex index: {:?}", min_idx);
        //}
        IdxEdge(min_idx, max_idx)
    }
}

impl From<Edge> for IdxEdge {
    fn from(edge: Edge) -> Self {
        IdxEdge::new(edge.0.idx, edge.1.idx)
    }
}

impl PartialEq for IdxEdge {
    fn eq(&self, other: &Self) -> bool {
        // WARN: This might match wrongly if the local indexing is not consistent.
        // Therfore, if any of the edges checked is formed from an intersection at its
        // turn, we can't ensure equality
        if matches!(self.0, PrimitiveIdx::Local(_))
            || matches!(self.1, PrimitiveIdx::Local(_))
            || matches!(other.0, PrimitiveIdx::Local(_))
            || matches!(other.1, PrimitiveIdx::Local(_))
        {
            panic!(
                "Cannot compare edges with local indexing, as the local indexing might not be consistent. Edge A: {:?}, Edge B: {:?}",
                self, other
            );
        }

        // NOTE: Due to the construction order of IdxEdge in IdxEdge::new(..),
        // (V0,V1) will only be encoded as (V0, V1), not (V1, V0)
        // Instead of (self.0 == other.0 && self.1 == other.1) || (self.0 == other.1 && self.1 == other.0), check:
        self.0 == other.0 && self.1 == other.1
    }
}
impl Eq for IdxEdge {}

#[derive(Debug, Clone, Copy, Hash)]
pub struct IdxIntersection(pub IdxEdge, pub IdxEdge);

impl IdxIntersection {
    pub fn new(idx_edge_a: IdxEdge, idx_edge_b: IdxEdge) -> Self {
        match idx_edge_a.0.cmp(&idx_edge_b.0) {
            std::cmp::Ordering::Less => IdxIntersection(idx_edge_a, idx_edge_b),
            std::cmp::Ordering::Greater => IdxIntersection(idx_edge_b, idx_edge_a),
            std::cmp::Ordering::Equal => {
                match idx_edge_a.1.cmp(&idx_edge_b.1) {
                    std::cmp::Ordering::Less => IdxIntersection(idx_edge_a, idx_edge_b),
                    std::cmp::Ordering::Greater => IdxIntersection(idx_edge_b, idx_edge_a),
                    std::cmp::Ordering::Equal => {
                        panic!(
                            "Edges are the same, cannot create intersection: {:?} and {:?}",
                            idx_edge_a, idx_edge_b
                        );
                        //IdxIntersection(idx_edge_a, idx_edge_b)
                    } // Edges are the same, order doesn't matter
                }
            }
        }
    }
}

impl PartialEq for IdxIntersection {
    fn eq(&self, other: &Self) -> bool {
        let IdxIntersection(idx_edge_a1, idx_edge_b1) = self;
        let IdxIntersection(idx_edge_a2, idx_edge_b2) = other;
        // NOTE: Due to the construction order of IdxIntersection and IdxEdge,
        // there will only be one intersection encoding for
        // intersection of Edge(U0,U1) and Edge(V0, V1)
        idx_edge_a1 == idx_edge_a2 && idx_edge_b1 == idx_edge_b2
    }
}
impl Eq for IdxIntersection {}

#[derive(Debug, Clone, Copy, Hash)]
pub enum PrimitiveIdx {
    Global(usize),
    Local(usize),
}

impl Default for PrimitiveIdx {
    fn default() -> Self {
        PrimitiveIdx::Local(usize::MAX)
    }
}

impl PartialEq for PrimitiveIdx {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (PrimitiveIdx::Global(idx_self), PrimitiveIdx::Global(idx_other)) => {
                idx_self == idx_other
            }
            (PrimitiveIdx::Local(idx_self), PrimitiveIdx::Local(idx_other)) => {
                idx_self == idx_other
            }
            _ => false,
        }
    }
}
impl Eq for PrimitiveIdx {}

impl Ord for PrimitiveIdx {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (PrimitiveIdx::Global(idx_self), PrimitiveIdx::Global(idx_other)) => {
                idx_self.cmp(idx_other)
            }
            (PrimitiveIdx::Local(idx_self), PrimitiveIdx::Local(idx_other)) => {
                idx_self.cmp(idx_other)
            }
            (PrimitiveIdx::Local { .. }, PrimitiveIdx::Global(_)) => std::cmp::Ordering::Less,
            (PrimitiveIdx::Global(_), PrimitiveIdx::Local { .. }) => std::cmp::Ordering::Greater,
        }
    }
}
impl PartialOrd for PrimitiveIdx {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Deref for PrimitiveIdx {
    type Target = usize;
    fn deref(&self) -> &Self::Target {
        match self {
            PrimitiveIdx::Local(idx) => idx,
            PrimitiveIdx::Global(idx) => idx,
        }
    }
}

#[derive(Clone, Debug, Copy, Default)]
pub struct Vertex {
    pub value: Point3<f64>,
    pub idx: PrimitiveIdx,
    // If vertex is produced from the intersection of two edges, register them here
    // NOTE: PrimitiveIdx::Local <=> from (IdxEdge, IdxEdge)
    // => Should perhaps encapsulate this information inside the PrimitiveIdx instead to have the
    // data hermetic and resolved by the compiler, instead of matching patterns
    pub from: Option<IdxIntersection>,
}

impl Hash for Vertex {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.idx.hash(state);
    }
}

// NOTE: Assume the vertex indexing is allocated properly such that
// no 2 distinct vertices will have the same idx
// WARN: That doesn't seem to be the case in the Wavefront file =>
// Must run collision detection and avoid allocating distinct PrimitiveIdx for the same vertex
// position. Not so critical for the norms
// WARN: When reallocating indexing, need to have clear barier in the code where the old and new
// local indexing is used => Assume any matching between PrimitiveIdx::Local can't be done
// reliably, hence the intersection is checked instead since `Local index |-> Intersection`
impl PartialEq for Vertex {
    fn eq(&self, other: &Self) -> bool {
        if self.idx == other.idx {
            true
        }
        //else if matches!(self.idx, PrimitiveIdx::Global(_)) || matches!(other.idx, PrimitiveIdx::Global(_)) {
        //    // If indexing didn't match, but one of the indexes are resolved globally, we don't
        //    // proceed checking intersections
        //    false
        //}
        else if let Some(idx_inter_a) = self.from
            && let Some(idx_inter_b) = other.from
        {
            idx_inter_a == idx_inter_b
        } else {
            false
        }
    }
}
impl Eq for Vertex {}

impl From<Vertex> for Point3<f64> {
    fn from(vertex: Vertex) -> Self {
        vertex.value
    }
}

impl Vertex {
    pub fn new(idx: usize, verts: &Vec<Point3<f64>>) -> Self {
        Vertex {
            value: verts[idx],
            idx: PrimitiveIdx::Global(idx),
            from: None,
        }
    }

    pub fn from_edges(edges: &Vec<Edge>) -> Vec<Self> {
        let mut verts = Vec::new();
        for edge in edges {
            if !verts.contains(&edge.0) {
                verts.push(edge.0);
            }
            if !verts.contains(&edge.1) {
                verts.push(edge.1);
            }
        }
        verts
    }
}

#[derive(Clone, Debug, Copy)]
pub struct Normal {
    pub value: Dir3,
    pub idx: PrimitiveIdx,
}

impl Normal {
    pub fn new(idx: usize, norms: &Vec<Dir3>) -> Self {
        Normal {
            value: norms[idx],
            idx: PrimitiveIdx::Global(idx),
        }
    }
}

pub struct Polygon {
    pub verts: Vec<Vertex>,
}

impl Polygon {
    pub fn edges(&self) -> Vec<Edge> {
        let mut edges = Vec::new();
        for i in 0..self.verts.len() {
            let v1 = self.verts[i].clone();
            let v2 = self.verts[(i + 1) % self.verts.len()].clone();
            edges.push(Edge(v1, v2));
        }
        edges
    }

    pub fn triangulate(&self) -> Vec<IdxTriangle> {
        let mut triangles = Vec::new();
        let mut idx = 0;
        for i in 1..(self.verts.len() - 1) {
            triangles.push(IdxTriangle::new(
                vec![self.verts[0], self.verts[i], self.verts[i + 1]],
                None,
                PrimitiveIdx::Local(idx),
            ));
            idx += 1;
        }
        triangles
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Point3;

    use super::*;

    #[test]
    fn test_edge_equality() {
        let v1 = Vertex {
            value: Point3::new(0.0, 0.0, 0.0),
            idx: PrimitiveIdx::Global(0),
            from: None,
        };
        let v2 = Vertex {
            value: Point3::new(1.0, 0.0, 0.0),
            idx: PrimitiveIdx::Global(1),
            from: None,
        };
        let edge1 = Edge(v1, v2);
        let edge2 = Edge(v2, v1);
        assert_eq!(edge1, edge2);
    }

    #[test]
    fn test_vertex_equality() {
        // Identity based on index
        let v1 = Vertex {
            value: Point3::new(0.0, 0.0, 0.0),
            idx: PrimitiveIdx::Global(0),
            from: None,
        };
        let v2 = Vertex {
            value: Point3::new(0.0, 0.0, 0.0),
            idx: PrimitiveIdx::Global(0),
            from: None,
        };
        assert_eq!(v1, v2);

        // Identity based on index, disregarding actual position
        let v1 = Vertex {
            value: Point3::new(0.0, 0.0, 0.0),
            idx: PrimitiveIdx::Global(0),
            from: None,
        };
        let v2 = Vertex {
            value: Point3::new(1.0, 0.0, 0.0),
            idx: PrimitiveIdx::Global(0),
            from: None,
        };
        assert_eq!(v1, v2);

        // Identity based on local index, disregarding actual position
        let v1 = Vertex {
            value: Point3::new(0.0, 0.0, 0.0),
            idx: PrimitiveIdx::Local(0),
            from: None,
        };
        let v2 = Vertex {
            value: Point3::new(1.0, 0.0, 0.0),
            idx: PrimitiveIdx::Local(0),
            from: None,
        };
        assert_eq!(v1, v2);

        // Check equality based on intersection
        let v1 = Vertex {
            value: Point3::new(0.0, 0.0, 0.0),
            idx: PrimitiveIdx::Local(0),
            from: Some(IdxIntersection::new(
                IdxEdge::new(PrimitiveIdx::Global(10), PrimitiveIdx::Global(11)),
                IdxEdge::new(PrimitiveIdx::Global(20), PrimitiveIdx::Global(21)),
            )),
        };
        // Check based on intersection
        let v2 = Vertex {
            value: Point3::new(1.0, 0.0, 0.0),
            idx: PrimitiveIdx::Local(1),
            from: Some(IdxIntersection::new(
                IdxEdge::new(PrimitiveIdx::Global(10), PrimitiveIdx::Global(11)),
                IdxEdge::new(PrimitiveIdx::Global(20), PrimitiveIdx::Global(21)),
            )),
        };
        assert_eq!(v1, v2);

        // Not equal vertices
        let v1 = Vertex {
            value: Point3::new(0.0, 0.0, 0.0),
            idx: PrimitiveIdx::Local(0),
            from: Some(IdxIntersection::new(
                IdxEdge::new(PrimitiveIdx::Global(10), PrimitiveIdx::Global(11)),
                IdxEdge::new(PrimitiveIdx::Global(20), PrimitiveIdx::Global(21)),
            )),
        };
        // Check based on intersection
        let v2 = Vertex {
            value: Point3::new(1.0, 0.0, 0.0),
            idx: PrimitiveIdx::Local(1),
            from: Some(IdxIntersection::new(
                IdxEdge::new(PrimitiveIdx::Global(10), PrimitiveIdx::Global(11)),
                IdxEdge::new(PrimitiveIdx::Global(20), PrimitiveIdx::Global(22)),
            )),
        };
        assert_ne!(v1, v2);
    }

    #[test]
    #[should_panic]
    fn test_vertex_local_intersection_invalid() {
        // Not equal vertices
        let v1 = Vertex {
            value: Point3::new(0.0, 0.0, 0.0),
            idx: PrimitiveIdx::Local(0),
            from: Some(IdxIntersection::new(
                IdxEdge::new(PrimitiveIdx::Global(10), PrimitiveIdx::Global(11)),
                IdxEdge::new(PrimitiveIdx::Global(20), PrimitiveIdx::Global(21)),
            )),
        };
        // Check based on intersection
        let v2 = Vertex {
            value: Point3::new(1.0, 0.0, 0.0),
            idx: PrimitiveIdx::Local(1),
            from: Some(IdxIntersection::new(
                // NOTE: Intersection from locally resolved edge not permited
                IdxEdge::new(PrimitiveIdx::Local(10), PrimitiveIdx::Global(11)),
                IdxEdge::new(PrimitiveIdx::Global(20), PrimitiveIdx::Global(21)),
            )),
        };
        assert_ne!(v1, v2);
    }

    #[test]
    fn test_idx_triangle_new_and_edges() {
        let verts = vec![
            Vertex {
                value: Point3::new(0.0, 0.0, 0.0),
                idx: PrimitiveIdx::Global(0),
                from: None,
            },
            Vertex {
                value: Point3::new(1.0, 0.0, 0.0),
                idx: PrimitiveIdx::Global(1),
                from: None,
            },
            Vertex {
                value: Point3::new(0.0, 1.0, 0.0),
                idx: PrimitiveIdx::Global(2),
                from: None,
            },
        ];
        let tri = IdxTriangle::new(verts.clone(), None, PrimitiveIdx::Global(0));
        let edges = tri.edges();
        assert_eq!(edges.len(), 3);
        assert!(edges.contains(&Edge(verts[0], verts[1])));
        assert!(edges.contains(&Edge(verts[1], verts[2])));
        assert!(edges.contains(&Edge(verts[2], verts[0])));
    }
}
