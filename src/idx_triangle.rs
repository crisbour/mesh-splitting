use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
};

use anyhow::{Result, anyhow};
use assert_approx_eq::assert_approx_eq;
use log::{debug, info, trace, warn};
use nalgebra::{Point3, Unit, Vector3};

use crate::{
    Collide, Segment, Split, SplitEdges, Triangle,
    primitives::{Edge, IdxEdge, IdxIntersection, Normal, Polygon, PrimitiveIdx, Vertex},
};

type Dir3 = Unit<Vector3<f64>>;

#[derive(Debug, Clone, Default)]
pub struct IdxTriangle {
    pub idx: PrimitiveIdx,
    pub verts: [Vertex; 3],
    pub norms: Option<[Normal; 3]>,
    pub flat: bool,
}

impl IdxTriangle {
    // TODO: Maybe need the parameters to be [T; 3], such that compiler checks the lengths instead
    pub fn new(verts: Vec<Vertex>, norms: Option<Vec<Normal>>, idx: PrimitiveIdx) -> Self {
        assert!(
            verts.len() == 3,
            "Expected exactly 3 vertices to construct a triangle"
        );
        let mut flat = true;
        if let Some(ref idx_norms) = norms {
            assert!(
                idx_norms.len() == 3,
                "Expected exactly 3 normals to construct a smooth triangle"
            );
            let v0 = idx_norms[0].value.into_inner();
            let v1 = idx_norms[1].value.into_inner();
            let v2 = idx_norms[2].value.into_inner();

            if (v1-v0).abs().sum() < 1e-10 && (v2-v0).abs().sum() < 1e-10 {
                flat = true;
            }
        }
        IdxTriangle {
            idx,
            verts: [verts[0], verts[1], verts[2]],
            norms: norms.map(|n| [n[0], n[1], n[2]]),
            flat,
        }
    }

    pub fn into_global(mut self, global_idx: usize) -> Self {
        self.idx = PrimitiveIdx::Global(global_idx);
        self
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

    /// Get barycentric_coords of a point in the triangle.
    /// Derived from "Real-Time Collision Detection" - Chapter 3.4
    fn barycentric_coords(&self, point: &Point3<f64>) -> Result<(f64, f64, f64)> {
        const EPS: f64 = 1e-9;
        let verts = self.verts;
        let v0 = verts[1].value - verts[0].value;
        let v1 = verts[2].value - verts[0].value;
        let v2 = *point - verts[0].value;
        let d00 = v0.dot(&v0);
        let d01 = v0.dot(&v1);
        let d11 = v1.dot(&v1);
        let d20 = v2.dot(&v0);
        let d21 = v2.dot(&v1);
        let denom = (d00 * d11) - (d01 * d01);
        if denom.abs() < f64::EPSILON {
            // Degenerate triangle (colinear edges), can't calculate barycentric coordinates.
            return Err(anyhow!(
                "Trying to calculate barycentric coordinates for a degenerate triangle."
            ));
        }
        let mut v = ((d11 * d20) - (d01 * d21)) / denom;
        let mut w = ((d00 * d21) - (d01 * d20)) / denom;
        let mut u = 1.0 - v - w;
        if u < -EPS || u > 1.0 + EPS || v < -EPS || v > 1.0 + EPS || w < -EPS || w > 1.0 + EPS {
            // Invalid barycentric coordinates, point outside the triangle
            return Err(anyhow!("Trying to calculate barycentric coordinates ({},{},{}) for a point outside the triangle.", u,v,w));
        }
        // Clamp barycentric coordinates to [0,1] to handle numerical precision issues
        u = u.max(0.0).min(1.0);
        v = v.max(0.0).min(1.0);
        w = w.max(0.0).min(1.0);
        assert!(
            (u + v + w - 1.0).abs() < EPS,
            "Barycentric coordinates do not sum to 1: u={}, v={}, w={}", u, v, w
        );
        Ok((u, v, w))
    }

    fn dir_at(&self, bary_coords: (f64, f64, f64)) -> Result<Dir3> {
        let idx_norms = self.norms.as_ref().expect("Expected normals for smooth triangle");
        let norm_value: Vector3<f64> = idx_norms[0].value.as_ref() * bary_coords.0
            + idx_norms[1].value.as_ref() * bary_coords.1
            + idx_norms[2].value.as_ref() * bary_coords.2;
        // FIXME: Perhaps only do this in debug mode for performance reasons
        assert_approx_eq!(norm_value.norm(), 1.0);
        Ok(Dir3::new_normalize(norm_value))
    }


    pub fn from_edges(&self, edges: Vec<Edge>) -> Vec<IdxTriangle> {
        if edges.is_empty() {
            return Vec::new();
        }

        // Group edges to form triangles
        let mut split_edges = SplitEdges::new(vec![]);
        split_edges.edges.extend(edges);
        let mut tris: Vec<IdxTriangle> = split_edges.into();

        // Compute norms if they exist in the original triangle
        if let Some(ref idx_norms) = self.norms {
            if self.flat {
                // If the original triangle is flat, we can directly assign the same normal to all new triangles without interpolation
                for tri in tris.iter_mut() {
                    tri.norms = Some([idx_norms[0], idx_norms[1], idx_norms[2]]);
                }
                return tris;
            }

            let mut norm_idx = 0;
            let tris_norms: Vec<_> = tris.iter().map(|tri| {
                // Norms for all vertices in triangle `tri`
                let norms = tri.verts.map(|v| {
                    if self.verts.contains(&v) {
                        // If the vertex of the new triangle matches one of the original triangle vertices, we can directly use the corresponding normal without interpolation
                        let vert_idx = self.verts.iter().position(|&orig_v| orig_v == v).unwrap();
                        idx_norms[vert_idx]
                    } else {
                        // Populate norms based on barycentric coordinates interpolation
                        let bc = self.barycentric_coords(&v.value).unwrap();
                        let norm = Normal { value: self.dir_at(bc).unwrap(), idx: PrimitiveIdx::Local(norm_idx)};
                        norm_idx += 1;
                        norm
                    }
                });
                norms
            }).collect();

            // Assign norms to the triangles
            for (tri, norms) in tris.iter_mut().zip(tris_norms){
                tri.norms = Some(norms);
            };
        }

        tris
    }

    pub fn rename_norms(&mut self, norm_idx_map: &HashMap<PrimitiveIdx, PrimitiveIdx>) {
        if let Some(ref mut idx_norms) = self.norms {
            for norm in idx_norms.iter_mut() {
                if let Some(new_idx) = norm_idx_map.get(&norm.idx) {
                    norm.idx = *new_idx;
                }
            }
        }
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
        IdxTriangle::new(poly.verts, None, PrimitiveIdx::Local(usize::MAX))
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

        let split_convex_polygon = SplitEdges::new(other.intersect(&new_edges));

        info!("Intersection edges: {:?}", split_convex_polygon.edges.iter().map(|e| Into::<IdxEdge>::into(e.clone())).collect::<Vec<_>>());

        split_convex_polygon.triangulate()
    }

    #[inline]
    fn split(&self, other: SplitEdges) -> (Self::Inst, SplitEdges) {
        trace!("Splitting {:?} with: {:?}", self, other);
        let mut split_edges = other.edges.clone();
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
            for seg_edge in other.outer_edges().iter() {
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

        trace!("Colinear vertices: {:?}", colinear_blacklist);

        // By now split_edges must contain all vertices in some edges from inter_edges and newly
        // constructed edges
        let verts = Vertex::from_edges(&split_edges);

        // The `new_edges` already form the intersection convex polygon
        for vert in verts.iter() {
            for other_vert in verts.iter() {
                if vert == other_vert {
                    continue;
                }
                if colinear_blacklist.get(&vert.idx).map_or(false, |set| set.contains(&other_vert.idx)) {
                    continue;
                }
                if split_edges.contains(&Edge(*vert, *other_vert)) {
                    // Skip edge which inside the splitting set, but not on the outer edges
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
        for edge in other.outer_edges().iter() {
            if diff_verts.contains(&edge.0) && diff_verts.contains(&edge.1) {
                if colinear_blacklist.get(&edge.0.idx).map_or(false, |set| set.contains(&edge.1.idx)) {
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

        info!("Split triangles: {:?}", tris_diff_vec.iter().map(|tri| tri.verts.map(|v| *v.idx)).collect::<Vec<_>>());

        (tris_diff_vec, other)
    }
}

impl Collide<Vec<Edge>> for IdxTriangle {
    fn overlap(&self, edges: &Vec<Edge>) -> bool {
        edges.iter().any(|edge| self.overlap(edge))
    }
}

impl Split<Vec<Edge>, Vec<Edge>> for IdxTriangle {
    type Inst = Vec<Edge>;
    /// Extract the segments that are overlapping the target triangle
    fn intersect(&self, edges: &Vec<Edge>) -> Vec<Edge> {
        trace!("Splitting triangle with segments: {:?}", edges);

        // Convert to simple spatial triangle
        let tri = self.tri();

        let mut verts: Vec<Vertex> = Vec::new();
        edges.iter().for_each(|Edge(v0, v1)| {
            if !verts.contains(&v0) && tri.overlap(&v0.value) {
                verts.push(*v0);
            }
            if !verts.contains(v1) && tri.overlap(&v0.value) {
                verts.push(*v1);
            }
        });

        // Only include segments that have both ends overlapping the triangle
        let edges = edges
            .iter()
            .filter(|Edge(v0, v1)|
                tri.overlap(&v0.value) && tri.overlap(&v1.value)
            )
            .cloned()
            .collect();
        trace!("Segments after intersection with triangle: {:?}", edges);

        edges
    }

    fn split(&self, _edges: Vec<Edge>) -> (Vec<Edge>, Vec<Edge>) {
        panic!("Not implemented yet");
    }
}
