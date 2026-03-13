use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
};

use log::{debug, info, trace, warn};
use nalgebra::{Unit, Vector3};

use crate::{
    Collide, Segment, Split, SplitEdges, Triangle,
    primitives::{Edge, IdxEdge, IdxIntersection, Normal, Polygon, PrimitiveIdx, Vertex},
};

type Dir3 = Unit<Vector3<f64>>;

#[derive(Debug, Clone, Default)]
pub struct IdxTriangle {
    pub verts: [Vertex; 3],
    pub norms: Option<[Normal; 3]>,
}

impl IdxTriangle {
    // TODO: Maybe need the parameters to be [T; 3], such that compiler checks the lengths instead
    pub fn new(verts: Vec<Vertex>, norms: Option<Vec<Normal>>) -> Self {
        assert!(
            verts.len() == 3,
            "Expected exactly 3 vertices to construct a triangle"
        );
        if let Some(ref idx_norms) = norms {
            assert!(
                idx_norms.len() == 3,
                "Expected exactly 3 normals to construct a smooth triangle"
            );
        }
        IdxTriangle {
            verts: [verts[0], verts[1], verts[2]],
            norms: norms.map(|n| [n[0], n[1], n[2]]),
        }
    }

    pub fn edges(&self) -> Vec<Edge> {
        vec![
            Edge(self.verts[0].clone(), self.verts[1].clone()),
            Edge(self.verts[1].clone(), self.verts[2].clone()),
            Edge(self.verts[2].clone(), self.verts[0].clone()),
        ]
    }

    pub fn tri(&self) -> Triangle {
        self.clone().into()
    }

    pub fn from_edges(&self, edges: Vec<Edge>) -> Vec<IdxTriangle> {
        let mut edge_map: HashMap<IdxEdge, Edge> = HashMap::new();
        for edge in edges.iter() {
            edge_map.insert(edge.clone().into(), *edge);
        }
        let verts = Vertex::from_edges(&edges);
        assert!(
            self.verts.iter().all(|v| verts.contains(v)),
            "All triangle vertices must be present in the edges"
        );
        let mut triangles = Vec::new();
        for i in 0..verts.len() {
            for j in (i + 1)..verts.len() {
                for k in (j + 1)..verts.len() {
                    let v1 = verts[i];
                    let v2 = verts[j];
                    let v3 = verts[k];
                    if edge_map.contains_key(&IdxEdge::new(v1.idx, v2.idx))
                        && edge_map.contains_key(&IdxEdge::new(v2.idx, v3.idx))
                        && edge_map.contains_key(&IdxEdge::new(v3.idx, v1.idx))
                    {
                        triangles.push(IdxTriangle::new(
                            vec![v1, v2, v3],
                            None,
                            //self.norms.map(|n| vec![n[i], n[j], n[k]]),
                        ));
                    }
                }
            }
        }
        // TODO: Check that all edges are used
        triangles
    }
}

impl Into<Triangle> for IdxTriangle {
    fn into(self) -> Triangle {
        let tri_verts = self
            .verts
            .iter()
            .map(|vert| vert.value)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        if let Some(ref idx_norms) = self.norms {
            let plane_norm = Dir3::new_normalize(
                idx_norms
                    .iter()
                    .fold(Vector3::zeros(), |acc, &norm| acc + norm.value.into_inner()),
            );
            Triangle::new_with_norm(tri_verts, plane_norm)
        } else {
            Triangle::new(tri_verts)
        }
    }
}

impl From<Vec<Edge>> for Polygon {
    fn from(edges: Vec<Edge>) -> Self {
        let mut verts = vec![edges[0].0];
        let end_vertex = edges[0].0;
        let mut next_vertex = edges[0].1;

        let mut visited = BTreeSet::new();
        visited.insert((edges[0].0.idx, edges[0].1.idx));

        while next_vertex != end_vertex {
            verts.push(next_vertex);
            if let Some(edge) = edges
                .iter()
                .find(|e| !visited.contains(&(e.0.idx, e.1.idx)) && e.0 == next_vertex)
            {
                next_vertex = edge.1;
                visited.insert((edge.0.idx, edge.1.idx));
            } else if let Some(edge) = edges
                .iter()
                .find(|e| !visited.contains(&(e.0.idx, e.1.idx)) && e.1 == next_vertex)
            {
                next_vertex = edge.0;
                visited.insert((edge.0.idx, edge.1.idx));
            } else {
                panic!("Edges do not form a closed loop");
            }
        }

        // Check that all edges have been visited
        for edge in edges.iter() {
            assert!(
                visited.contains(&(edge.0.idx, edge.1.idx)),
                "Not all edges have been used in constructing the polygon.\nPolygon: {:?}\nEdge unused {:?}\nFrom edges: {:?}",
                verts.iter().map(|v| v.idx).collect::<Vec<_>>(),
                edge,
                edges.iter().map(|e| e.clone().into()).collect::<Vec<IdxEdge>>()
            );
        }

        Polygon { verts }
    }
}

impl From<IdxTriangle> for Polygon {
    fn from(tri: IdxTriangle) -> Self {
        Polygon {
            verts: tri.verts.to_vec(),
        }
    }
}

impl From<Polygon> for IdxTriangle {
    fn from(poly: Polygon) -> Self {
        assert!(
            poly.verts.len() == 3,
            "Expected exactly 3 vertices to construct a triangle"
        );
        IdxTriangle::new(poly.verts, None)
    }
}

impl Collide<IdxTriangle> for IdxTriangle {
    fn overlap(&self, other: &IdxTriangle) -> bool {
        let tri_self: Triangle = self.clone().into();
        let tri_other: Triangle = other.clone().into();
        tri_self.overlap(&tri_other)
    }
}

impl Collide<Edge> for IdxTriangle {
    // WARN: This assumes Edge and IdxTriangle are coplanar
    fn overlap(&self, edge: &Edge) -> bool {
        // Edge must intersect a triangle edge without being fuzzy,
        // or two intersections exists.
        // Two fuzzy intersections are not valid if they are very close to one another
        //   - Avoid edge passing through vertex to be considered intersecting
        let seg = Segment::new(edge.0.value, edge.1.value);
        let inters: Vec<(_, _)> = self
            .edges()
            .iter()
            .map(|tri_edge| Segment::new(tri_edge.0.value, tri_edge.1.value))
            .flat_map(|tri_seg| tri_seg.intersect_with_eps(&seg, 1e-9))
            .collect();

        // FIXME: Reduce intersection identified close to triangle vertex
        // since that can result in 3 intersections, since 2 of them represent the same point
        assert!(inters.len() <= 2, "Expected at most 2 intersections");

        // Check if any intersection is certain
        for (_inter, fuzzy) in inters.iter() {
            if !fuzzy {
                return true;
            }
        }

        // Check if two fuzzy intersections are sufficiently far apart to be considered valid
        if inters.len() == 2 {
            let (inter1, _) = inters[0];
            let (inter2, _) = inters[1];
            let diff = inter1 - inter2;
            if diff.abs().sum() > 1e-9 { true } else { false }
        } else {
            false
        }
    }
}

impl Split<IdxTriangle, SplitEdges> for IdxTriangle {
    // TODO: Also need to solve non-coplanar triangles that intersect forming one and only
    // intersecting segment
    type Inst = Vec<Self>;

    fn intersect(&self, other: &IdxTriangle) -> SplitEdges {
        trace!("Check intersections of {:?} with {:?}", self, other);

        let mut local_verts: Vec<Vertex> = vec![];
        let mut local_verts_valid: Vec<bool> = vec![];

        let mut new_edges: Vec<Edge> = Vec::new();

        let mut verts_in = vec![];
        let mut verts_out = vec![];
        let mut edges_v = vec![];

        // Resolve vertices and segments inside the triangle
        for &v in &other.verts {
            if self.tri().vertex_in(v.into()) {
                verts_in.push(v);
            } else {
                verts_out.push(v);
            }
        }

        // Add segments between vertices inside the triangle, as they will be part of the new triangles
        for i in 0..verts_in.len() {
            for j in (i + 1)..verts_in.len() {
                new_edges.push(Edge(verts_in[i], verts_in[j]));
            }
        }
        trace!("Edges between vertices inside the triangle: {:?}", new_edges.iter().map(|e| Into::<IdxEdge>::into(e.clone())).map(|IdxEdge(v0,v1)| (*v0, *v1)).collect::<Vec<_>>());

        trace!("Verts in: {:?}, Verts out: {:?}", verts_in, verts_out);
        // Get all segments that might intersect the triangle
        match verts_in.len() {
            3 => {
                assert_eq!(verts_out.len(), 0);
                // No vertices outside. The other triangle is fully contained in this one.
            }
            2 => {
                assert_eq!(verts_out.len(), 1);
                edges_v.push((Some(verts_in[0]), Edge(verts_out[0], verts_in[0])));
                edges_v.push((Some(verts_in[1]), Edge(verts_out[0], verts_in[1])));
            }
            1 => {
                assert_eq!(verts_out.len(), 2);
                edges_v.push((Some(verts_in[0]), Edge(verts_out[0], verts_in[0])));
                edges_v.push((Some(verts_in[0]), Edge(verts_out[1], verts_in[0])));
                edges_v.push((None, Edge(verts_out[0], verts_out[1])));
            }
            0 => {
                assert_eq!(verts_out.len(), 3);
                edges_v.push((None, Edge(verts_out[0], verts_out[1])));
                edges_v.push((None, Edge(verts_out[1], verts_out[2])));
                edges_v.push((None, Edge(verts_out[2], verts_out[0])));
            }
            _ => unreachable!(),
        }

        trace!("Segments to check intersections: {:?}", edges_v);

        // Find all intersections and segments that need to be included in the new triangles
        #[derive(Debug, Clone)]
        struct Intersection {
            value: Vertex,
            fuzzy: bool,
        }

        let mut inters: Vec<Intersection> = Vec::new();

        // Find all intersections for each edge of the other triangle
        for (vertex_v, edge_v) in edges_v {
            let seg_v = Segment::new(edge_v.0.value, edge_v.1.value);
            trace!(
                "Vertex {:?} is inside the triangle. Segment: {:?}",
                vertex_v, seg_v
            );
            let mut edge_v_inters = Vec::new();
            let intersections: Vec<_> = self
                .edges()
                .iter()
                .map(|Edge(u1, u2)| {
                    (
                        Segment::new(u1.value, u2.value).intersect_with_eps(&seg_v, 1e-9),
                        u1,
                        u2,
                    )
                })
                .filter(|(intersection, _u1, _u2)| intersection.is_some())
                .map(|(intersection, u1, u2)| (intersection.unwrap(), *u1, *u2))
                .collect();

            trace!("Intersections: {:?}", intersections);

            let mut inter_indices = Vec::new();
            for (intersection, u1, u2) in intersections {
                // Decide the index of the intersection vertex. It can be either of the edge
                // vertices if the intersection is close enough, a new vertex or reuse a previous
                // constructed vertex
                //   - The global vertex idx is used instead of one is identified
                //   - The position is overwritten with the identified global vertex to avoid
                // opennings in the mesh
                let fuzzy = intersection.1;
                let (idx, inter_pos) = if (u1.value - intersection.0).abs().sum() < 1e-9 {
                    (u1.idx, u1.value)
                } else if (u2.value - intersection.0).abs().sum() < 1e-9 {
                    (u2.idx, u2.value)
                } else if (edge_v.0.value - intersection.0).abs().sum() < 1e-9 {
                    (edge_v.0.idx, edge_v.0.value)
                } else if (edge_v.1.value - intersection.0).abs().sum() < 1e-9 {
                    (edge_v.1.idx, edge_v.1.value)
                } else {
                    // Allocate a new index or reuse a local one if positions matches withing some
                    // norm1 eps error
                    let mut idx = PrimitiveIdx::Local(local_verts.len());
                    let mut value = intersection.0;
                    for &v in local_verts.iter() {
                        if (v.value - intersection.0).abs().sum() < 1e-9 {
                            idx = v.idx;
                            value = v.value;
                            break;
                        }
                    }
                    (idx, value)
                };

                if !inter_indices.contains(&idx) {
                    inter_indices.push(idx);
                    let inter_vertex = Vertex {
                        idx,
                        value: inter_pos,
                        // NOTE: Event though idx might be PrimitiveIdx::Global, we keep track of
                        // the fact this vertex is being used as a result of an intersection
                        // necessity
                        from: Some(IdxIntersection(Edge(u1, u2).into(), edge_v.into())),
                    };
                    if idx != u1.idx && idx != u2.idx && idx != edge_v.0.idx && idx != edge_v.1.idx
                    {
                        local_verts.push(inter_vertex);
                        local_verts_valid.push(false);
                    }
                    edge_v_inters.push(Intersection {
                        value: inter_vertex,
                        fuzzy,
                    });
                }
            }

            assert!(
                edge_v_inters.len() <= 2,
                "Expected at most 2 intersections, but {} where found.",
                edge_v_inters.len()
            );

            let inters_count = edge_v_inters.len();
            let edge_v_inters_validated: Vec<Intersection> = edge_v_inters
                .into_iter()
                .filter(|inter| vertex_v.is_some() || !inter.fuzzy || inters_count == 2)
                .collect();

            edge_v_inters_validated
                .iter()
                .for_each(|inter| match inter.value.idx {
                    PrimitiveIdx::Global(idx) => {
                        trace!("Intersection at global vertex index {}", idx);
                    }
                    PrimitiveIdx::Local(idx) => {
                        trace!("Intersection at local vertex index {}", idx);
                        local_verts_valid[idx] = true;
                    }
                });

            trace!("Intersections: {:?}", edge_v_inters_validated);

            trace!(
                "Vertices {:?} with validity {:?}",
                local_verts, local_verts_valid
            );

            if let Some(vertex_v) = vertex_v {
                trace!(
                    "Vertex {:?} is inside the triangle {:?}. Intersections: {}",
                    vertex_v,
                    self.verts.iter().map(|v| v.value),
                    edge_v_inters_validated.len()
                );
                assert_eq!(edge_v_inters_validated.len(), 1);
                new_edges.push(Edge(vertex_v, edge_v_inters_validated[0].value));
            } else {
                if edge_v_inters_validated.len() == 2 {
                    new_edges.push(Edge(
                        edge_v_inters_validated[0].value,
                        edge_v_inters_validated[1].value,
                    ));
                } else {
                    if edge_v_inters_validated.len() == 1 {
                        debug!("WARN: Outer intersection with eps error must be checked");
                    }
                }
            }
            inters.extend(edge_v_inters_validated);
        }

        // Deduplicate inter_idx among segs_inter
        let mut inter_indices = Vec::new();
        let mut i = 0;
        while i < inters.len() {
            if inter_indices.contains(&inters[i].value) {
                inters.remove(i);
            } else {
                inter_indices.push(inters[i].value);
                i += 1;
            }
        }
        trace!("Intersections: {:?}", inters);

        trace!("New segments after intersection: {:?}", new_edges);

        // Solving colinear vertices and generate segs that don't overlap each other
        // -------------------------------------------------------------------------
        let mut colinear_blacklist: BTreeMap<PrimitiveIdx, Vec<PrimitiveIdx>> = BTreeMap::new();
        // Populate edge entry with its extreme points
        for edge in self.edges() {
            let seg = Segment::new(edge.0.value, edge.1.value);
            let mut colinear_verts = Vec::new();
            colinear_verts.push(edge.0);
            colinear_verts.push(edge.1);
            for inter_vertex in inters.as_slice() {
                if seg.overlap(&inter_vertex.value.value) {
                    colinear_verts.push(inter_vertex.value);
                }
            }
            colinear_verts.sort_by(|v_idx_0, v_idx_1| {
                let alpha_0 = seg.parametric_dist(&v_idx_0.value);
                let alpha_1 = seg.parametric_dist(&v_idx_1.value);
                alpha_0
                    .partial_cmp(&alpha_1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }); // Sort by alpha
            for i in 0..colinear_verts.len() {
                let v_i = colinear_verts[i];
                if i > 0 {
                    let v_left = colinear_verts[i - 1];
                    if !new_edges.contains(&Edge(v_left, v_i)) {
                        new_edges.push(Edge(v_left, v_i));
                        trace!(
                            "Adding segment between colinear vertices: {:?} and {:?}",
                            v_left.idx, v_i.idx
                        );
                    }
                }
                if i < colinear_verts.len() - 1 {
                    let v_right = colinear_verts[i + 1];
                    if !new_edges.contains(&Edge(v_i, v_right)) {
                        new_edges.push(Edge(v_i, v_right));
                        trace!(
                            "Adding segment between colinear vertices: {:?} and {:?}",
                            v_right.idx, v_i.idx
                        );
                    }
                }
                for j in 0..colinear_verts.len() {
                    if ((j as isize) < (i as isize - 1)) || (j > i + 1) {
                        let v_j = colinear_verts[j];
                        let blacklist_j = colinear_blacklist.entry(v_j.idx).or_default();
                        if !blacklist_j.contains(&v_i.idx) {
                            blacklist_j.push(v_i.idx);
                        }
                        let blacklist_i = colinear_blacklist.entry(v_i.idx).or_default();
                        if !blacklist_i.contains(&v_j.idx) {
                            blacklist_i.push(v_j.idx);
                        }
                    }
                }
            }
        }

        trace!("New segments after handling colinearity: {:?}", new_edges.iter().map(|e| Into::<IdxEdge>::into(e.clone())).collect::<Vec<_>>());
        trace!("Colinear vertices: {:?}", colinear_blacklist);

        // TODO: Instead, we can complete the `new_segs` that would reprsent convex polygon
        let split_convex_polygon = SplitEdges(new_edges).intersect(&other);

        info!("Intersection edges: {:?}", split_convex_polygon.0.iter().map(|e| Into::<IdxEdge>::into(e.clone())).collect::<Vec<_>>());

        let convex_polygon = Polygon::from(split_convex_polygon.0.clone());

        let triangles = convex_polygon.triangulate();
        triangles.into()
    }

    #[inline]
    fn split(&self, other: SplitEdges) -> (Self::Inst, SplitEdges) {
        let inter_edges = other.0.clone();
        let mut split_edges = other.0.clone();
        let mut diff_edges = Vec::new();

        // Solving colinear vertices and generate segs that don't overlap each other
        // -------------------------------------------------------------------------
        let mut colinear_blacklist: BTreeMap<PrimitiveIdx, Vec<PrimitiveIdx>> = BTreeMap::new();
        // Populate edge entry with its extreme points
        for edge in self.edges() {
            let seg = Segment::new(edge.0.value, edge.1.value);
            let mut colinear_verts = Vec::new();
            colinear_verts.push(edge.0);
            colinear_verts.push(edge.1);
            for seg_edge in other.0.iter() {
                for vi in [seg_edge.0, seg_edge.1] {
                    // WARN: What about intersection which matches the original vertex of
                    // V triangle?
                    if let Some(IdxIntersection(u_edge, v_edge)) = vi.from {
                        if IdxEdge::from(edge) == u_edge || IdxEdge::from(edge) == v_edge {
                            colinear_verts.push(vi);
                        }
                    }
                }
            }
            colinear_verts.sort_by(|v_idx_0, v_idx_1| {
                let alpha_0 = seg.parametric_dist(&v_idx_0.value);
                let alpha_1 = seg.parametric_dist(&v_idx_1.value);
                alpha_0
                    .partial_cmp(&alpha_1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }); // Sort by alpha
            for i in 0..colinear_verts.len() {
                let v_i = colinear_verts[i];
                if i > 0 {
                    let v_left = colinear_verts[i - 1];
                    if !split_edges.contains(&Edge(v_left, v_i)) {
                        split_edges.push(Edge(v_left, v_i));
                        diff_edges.push(Edge(v_left, v_i));
                    }
                }
                if i < colinear_verts.len() - 1 {
                    let v_right = colinear_verts[i + 1];
                    if !split_edges.contains(&Edge(v_i, v_right)) {
                        split_edges.push(Edge(v_i, v_right));
                        diff_edges.push(Edge(v_i, v_right));
                    }
                }
                for j in 0..colinear_verts.len() {
                    if ((j as isize) < (i as isize - 1)) || (j > i + 1) {
                        let v_j = colinear_verts[j];
                        let blacklist_j = colinear_blacklist.entry(v_j.idx).or_default();
                        if !blacklist_j.contains(&v_i.idx) {
                            blacklist_j.push(v_i.idx);
                        }
                        let blacklist_i = colinear_blacklist.entry(v_i.idx).or_default();
                        if !blacklist_i.contains(&v_j.idx) {
                            blacklist_i.push(v_j.idx);
                        }
                    }
                }
            }
        }

        // By now split_edges must contain all vertices in some edges from inter_edges and newly
        // constructed edges
        let verts = Vertex::from_edges(&split_edges);

        // The `new_edges` already form the intersection convex polygon
        for vert in verts.iter() {
            for other_vert in verts.iter() {
                if vert == other_vert {
                    continue;
                }
                if colinear_blacklist.contains_key(&vert.idx) && colinear_blacklist[&vert.idx].contains(&other_vert.idx) {
                    continue;
                }
                let new_seg = Segment::new(vert.value, other_vert.value);
                let mut clash = false;
                for Edge(v0, v1) in split_edges.as_slice() {
                    match Segment::new(v0.value, v1.value).intersect_with_eps(&new_seg, 1e-9) {
                        Some((_intersect, false)) => {
                            clash = true;
                            break;
                        }
                        _ => {}
                    }
                }
                if !clash {
                    split_edges.push(Edge(*vert, *other_vert));
                    diff_edges.push(Edge(*vert, *other_vert));
                }
            }
        }

        // Include boundary edges from intersection edges
        let diff_verts = Vertex::from_edges(&diff_edges);
        for edge in inter_edges.iter() {
            if diff_verts.contains(&edge.0) && diff_verts.contains(&edge.1) {
                if colinear_blacklist.contains_key(&edge.0.idx) && colinear_blacklist[&edge.0.idx].contains(&edge.1.idx) {
                    // Skip edges that are blacklisted
                    continue;
                }
                if !diff_edges.contains(edge) {
                    split_edges.push(*edge);
                    diff_edges.push(*edge);
                } else {
                    warn!(
                        "Had not expect to find an edge from intersection already included in the difference set"
                    );
                }
            }
        }

        let tris_diff_vec = self.from_edges(diff_edges.clone());

        (tris_diff_vec, SplitEdges(inter_edges))
    }
}
