use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fmt::Display,
};

use anyhow::{Context, Result, anyhow};
use assert_approx_eq::assert_approx_eq;
use log::{info, trace, warn};
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

impl PartialEq for IdxTriangle {
    fn eq(&self, other: &Self) -> bool {
        let self_idxs = self.verts.map(|v| v.idx);
        let other_idxs = other.verts.map(|v| v.idx);
        self_idxs.iter().all(|v| other_idxs.contains(v))
            && other_idxs.iter().all(|v| self_idxs.contains(v))
    }
}

impl Display for IdxTriangle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "IdxTriangle {{ idx: {:?}, verts: [{}:{}, {}:{}, {}:{}] }}",
            self.idx,
            *self.verts[0].idx,
            self.verts[0].value,
            *self.verts[1].idx,
            self.verts[1].value,
            *self.verts[2].idx,
            self.verts[2].value,
        )
    }
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
            if idx_norms[0].idx == idx_norms[1].idx && idx_norms[1].idx == idx_norms[2].idx {
                flat = true;
            }
        }
        let a = verts[0].value - verts[1].value;
        let b = verts[1].value - verts[2].value;
        assert!(
            a.cross(&b).abs().sum() > f64::EPSILON,
            "Edges of triangle can't be colinear: {:?}",
            verts
        );
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

    pub fn edges(&self) -> [Edge; 3] {
        [
            Edge::new(self.verts[0].clone(), self.verts[1].clone()),
            Edge::new(self.verts[1].clone(), self.verts[2].clone()),
            Edge::new(self.verts[2].clone(), self.verts[0].clone()),
        ]
    }

    pub fn segs(&self) -> [Segment; 3] {
        self.edges().map(|e| e.into())
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
            return Err(anyhow!(
                "Trying to calculate barycentric coordinates ({},{},{}) for a point outside the triangle.",
                u,
                v,
                w
            ));
        }
        // Clamp barycentric coordinates to [0,1] to handle numerical precision issues
        u = u.max(0.0).min(1.0);
        v = v.max(0.0).min(1.0);
        w = w.max(0.0).min(1.0);
        assert!(
            (u + v + w - 1.0).abs() < EPS,
            "Barycentric coordinates do not sum to 1: u={}, v={}, w={}",
            u,
            v,
            w
        );
        Ok((u, v, w))
    }

    fn dir_at(&self, bary_coords: (f64, f64, f64)) -> Result<Dir3> {
        let idx_norms = self
            .norms
            .as_ref()
            .context("Expected normals for smooth triangle")?;
        let norm_value: Vector3<f64> = idx_norms[0].value.as_ref() * bary_coords.0
            + idx_norms[1].value.as_ref() * bary_coords.1
            + idx_norms[2].value.as_ref() * bary_coords.2;
        // FIXME: Perhaps only do this in debug mode for performance reasons
        assert_approx_eq!(norm_value.norm(), 1.0);
        Ok(Dir3::new_normalize(norm_value))
    }

    pub fn from_edges(&self, edges: Vec<Edge>) -> Result<Vec<IdxTriangle>> {
        if edges.is_empty() {
            return Ok(Vec::new());
        }

        // Group edges to form triangles
        let mut split_edges = SplitEdges::new(vec![]);
        split_edges.edges.extend(edges);
        let mut tris: Vec<IdxTriangle> = split_edges.try_into()?;

        // Compute norms if they exist in the original triangle
        if let Some(ref idx_norms) = self.norms {
            if self.flat {
                // If the original triangle is flat, we can directly assign the same normal to all new triangles without interpolation
                for tri in tris.iter_mut() {
                    tri.norms = Some([idx_norms[0], idx_norms[1], idx_norms[2]]);
                }
                return Ok(tris);
            }

            let mut norm_idx = 0;
            let tris_norms: Vec<_> = tris
                .iter()
                .map(|tri| {
                    // Norms for all vertices in triangle `tri`
                    let norms = tri.verts.map(|v| {
                        if self.verts.contains(&v) {
                            // If the vertex of the new triangle matches one of the original triangle vertices, we can directly use the corresponding normal without interpolation
                            let vert_idx =
                                self.verts.iter().position(|&orig_v| orig_v == v).unwrap();
                            idx_norms[vert_idx]
                        } else {
                            // Populate norms based on barycentric coordinates interpolation
                            let bc = self.barycentric_coords(&v.value).unwrap();
                            let norm = Normal {
                                value: self.dir_at(bc).unwrap(),
                                idx: PrimitiveIdx::Local(norm_idx),
                            };
                            norm_idx += 1;
                            norm
                        }
                    });
                    norms
                })
                .collect();

            // Assign norms to the triangles
            for (tri, norms) in tris.iter_mut().zip(tris_norms) {
                tri.norms = Some(norms);
            }
        }

        Ok(tris)
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

impl TryFrom<Vec<Edge>> for Polygon {
    type Error = anyhow::Error;
    fn try_from(edges: Vec<Edge>) -> Result<Self> {
        if edges.is_empty() {
            return Ok(Polygon { verts: vec![] });
        }
        if edges.len() < 3 {
            return Err(anyhow!(
                "Expected at least 3 edges to construct a polygon, but {} were provided.",
                edges.len()
            ));
        }

        // TODO: Add check that the degree of each vert is 2
        let mut verts = vec![edges[0].0];
        let end_vertex = edges[0].0;
        let mut next_vertex = edges[0].1;

        let mut visited: BTreeSet<IdxEdge> = BTreeSet::new();
        visited
            .insert(edges[0].try_into().context("Mark first edge used")?);

        while next_vertex != end_vertex {
            verts.push(next_vertex);
            let mut found = false;
            for edge in edges.iter() {
                let idx_edge: IdxEdge = edge
                    .clone()
                    .try_into()
                    .context(format!("Check next edge to be used: {:?}", edge))?;
                if !visited.contains(&idx_edge) && (edge.0 == next_vertex || edge.1 == next_vertex)
                {
                    found = true;
                    visited.insert(idx_edge);
                    if edge.0 == next_vertex {
                        next_vertex = edge.1;
                    } else {
                        next_vertex = edge.0;
                    }
                    break;
                }
            }
            if !found {
                return Err(anyhow!("Edges do not form a closed loop"));
            }
        }

        // Check that all edges have been visited
        for edge in edges.iter() {
            let idx_edge = IdxEdge::new(edge.0.idx, edge.1.idx)?;
            if !visited.contains(&idx_edge) {
                return Err(anyhow!(
                    "Not all edges have been used in constructing the polygon.\nPolygon: {:?}\nEdge unused {:?}\nFrom edges: {:?}",
                    verts.iter().map(|v| v.idx).collect::<Vec<_>>(),
                    idx_edge,
                    edges
                        .iter()
                        .map(|e| e.clone().try_into().unwrap())
                        .collect::<Vec<IdxEdge>>()
                ));
            }
        }

        Ok(Polygon { verts })
    }
}

impl From<IdxTriangle> for Polygon {
    fn from(tri: IdxTriangle) -> Self {
        Polygon {
            verts: tri.verts.to_vec(),
        }
    }
}

impl TryFrom<Polygon> for IdxTriangle {
    type Error = anyhow::Error;
    fn try_from(poly: Polygon) -> Result<Self> {
        if poly.verts.len() != 3 {
            Err(anyhow!(
                "Expected exactly 3 vertices to construct a triangle, but polygon has {}", poly.verts.len()
            ))
        } else {
            Ok(IdxTriangle::new(
                poly.verts,
                None,
                PrimitiveIdx::Local(usize::MAX),
            ))
        }
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

    fn intersect(&self, other: &IdxTriangle) -> Result<SplitEdges> {
        trace!("Check intersections of {} with {}", self, other);

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

        trace!("Verts in: {:?}, Verts out: {:?}", verts_in, verts_out);

        // Add segments between vertices inside the triangle, as they will be part of the new triangles
        for i in 0..verts_in.len() {
            for j in (i + 1)..verts_in.len() {
                new_edges.push(Edge::new(verts_in[i], verts_in[j]));
            }
        }
        trace!(
            "Edges between vertices inside the triangle: {:?}",
            new_edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()).unwrap())
                .map(|IdxEdge(v0, v1)| (*v0, *v1))
                .collect::<Vec<_>>()
        );

        // Get all segments that might intersect the triangle
        match verts_in.len() {
            3 => {
                assert_eq!(verts_out.len(), 0);
                // No vertices outside. The other triangle is fully contained in this one.
            }
            2 => {
                assert_eq!(verts_out.len(), 1);
                edges_v.push((Some(verts_in[0]), Edge::new(verts_out[0], verts_in[0])));
                edges_v.push((Some(verts_in[1]), Edge::new(verts_out[0], verts_in[1])));
            }
            1 => {
                assert_eq!(verts_out.len(), 2);
                edges_v.push((Some(verts_in[0]), Edge::new(verts_out[0], verts_in[0])));
                edges_v.push((Some(verts_in[0]), Edge::new(verts_out[1], verts_in[0])));
                edges_v.push((None, Edge::new(verts_out[0], verts_out[1])));
            }
            0 => {
                assert_eq!(verts_out.len(), 3);
                edges_v.push((None, Edge::new(verts_out[0], verts_out[1])));
                edges_v.push((None, Edge::new(verts_out[1], verts_out[2])));
                edges_v.push((None, Edge::new(verts_out[2], verts_out[0])));
            }
            _ => unreachable!(),
        }

        trace!(
            "Segments to check intersections: {:?}",
            edges_v
                .iter()
                .map(|(sv, e)| format!(
                    "(Vert: {:?}, Edge: {:?}), ",
                    sv.map(|v| v.idx),
                    IdxEdge::try_from(e.clone())
                ))
                .collect::<String>()
        );

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
                        from: Some(IdxIntersection(Edge::new(u1, u2).try_into().unwrap(), edge_v.try_into().unwrap())),
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

            edge_v_inters_validated.iter().for_each(|inter| {
                trace!("Intersection at {:?}", inter.value.idx);
                if let PrimitiveIdx::Local(idx) = inter.value.idx {
                    local_verts_valid[idx] = true;
                }
            });

            trace!("Validated Intersections: {:?}", edge_v_inters_validated);

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
                new_edges.push(Edge::new(vertex_v, edge_v_inters_validated[0].value));
            } else {
                if edge_v_inters_validated.len() == 2 {
                    new_edges.push(Edge::new(
                        edge_v_inters_validated[0].value,
                        edge_v_inters_validated[1].value,
                    ));
                } else {
                    if edge_v_inters_validated.len() == 1 {
                        warn!("Outer intersection with eps error must be checked");
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

        trace!(
            "New segments after intersection: {:?}",
            new_edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Vec<_>>()
        );

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
                // Ignore add vertex if it matches the target triangle (self) vertices
                if inter_vertex.value.idx == edge.0.idx {
                    // Replace original vertex, in order to preserve `from` field
                    colinear_verts[0] = inter_vertex.value;
                }
                else if inter_vertex.value.idx == edge.1.idx {
                    // Replace original vertex, in order to preserve `from` field
                    colinear_verts[1] = inter_vertex.value;
                }
                else {
                    if seg.overlap(&inter_vertex.value.value) {
                        colinear_verts.push(inter_vertex.value);
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
            trace!("Colinear vertices for edge {:?}: {:?}",
                IdxEdge::try_from(edge.clone())?,
                colinear_verts
                    .iter()
                    .map(|v| v)
                    .collect::<Vec<_>>()
            );
            for i in 0..colinear_verts.len() {
                let v_i = colinear_verts[i];
                if i > 0 {
                    let v_left = colinear_verts[i - 1];
                    if !new_edges.contains(&Edge::new(v_left, v_i)) {
                        new_edges.push(Edge::new(v_left, v_i));
                        trace!(
                            "Adding segment between colinear vertices: {:?} and {:?}",
                            v_left.idx, v_i.idx
                        );
                    }
                }
                if i < colinear_verts.len() - 1 {
                    let v_right = colinear_verts[i + 1];
                    if !new_edges.contains(&Edge::new(v_i, v_right)) {
                        new_edges.push(Edge::new(v_i, v_right));
                        trace!(
                            "Adding segment between colinear vertices: {:?} and {:?}",
                            v_right.idx, v_i.idx
                        );
                    }
                }
                // FIXME: This is not necessary for the intersection
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

        trace!(
            "New segments after handling colinearity: {:?}",
            new_edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Vec<_>>()
        );
        trace!("Colinear vertices: {:?}", colinear_blacklist);

        let split_convex_polygon = SplitEdges::new(other.intersect(&new_edges)?)
            .triangulate()?;

        info!(
            "Intersection edges: {:?} with outer: {:?}",
            split_convex_polygon
                .edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Vec<_>>(),
            split_convex_polygon
                .outer_edges()
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Vec<_>>()
        );

        Ok(split_convex_polygon)
    }

    #[inline]
    fn split(&self, other: SplitEdges) -> Result<(Self::Inst, SplitEdges)> {
        let mut split_edges = other.edges.clone();
        trace!(
            "Splitting {:?} with: {:?}",
            self,
            split_edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Result<Vec<_>>>()?
        );
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
            // TODO: Replace iteration through outer edges to interation through outer verts
            //let outer_verts = Vertex::from(other.outer_edges());
            for seg_edge in other.outer_edges().iter() {
                for vi in [seg_edge.0, seg_edge.1] {
                    // WARN: What about intersection which matches the original vertex of
                    // V triangle?
                    trace!("Checking if vertex {:?} is colinear with triangle edge {:?}",
                        vi, IdxEdge::try_from(edge.clone()).unwrap());
                    if let Some(IdxIntersection(u_edge, v_edge)) = vi.from {
                        if IdxEdge::try_from(edge).unwrap() == u_edge || IdxEdge::try_from(edge).unwrap() == v_edge {
                            // Prevent triangle original vertex to be double counted
                            if !colinear_verts.contains(&vi) {
                                colinear_verts.push(vi);
                            }
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
            trace!(
                "Colinear vertices for edge {:?}: {:?}",
                IdxEdge::try_from(edge.clone())?,
                colinear_verts
                    .iter()
                    .map(|v| v)
                    .collect::<Vec<_>>()
            );
            for i in 0..colinear_verts.len() {
                let v_i = colinear_verts[i];
                if i > 0 {
                    let v_left = colinear_verts[i - 1];
                    let new_edge = Edge::new(v_left, v_i);
                    if !split_edges.contains(&new_edge) {
                        split_edges.push(new_edge);
                        diff_edges.push(new_edge);
                    }
                }
                if i < colinear_verts.len() - 1 {
                    let v_right = colinear_verts[i + 1];
                    let new_edge = Edge::new(v_i, v_right);
                    if !split_edges.contains(&new_edge) {
                        split_edges.push(new_edge);
                        diff_edges.push(new_edge);
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
        trace!(
            "Diff edges after handling colinearity: {:?}",
            diff_edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Result<Vec<_>>>()?
        );
        trace!(
            "Split edges after handling colinearity: {:?}",
           split_edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Result<Vec<_>>>()?
        );

        trace!("Colinear vertices: {:?}", colinear_blacklist);

        // By now split_edges must contain all vertices in some edges from inter_edges and newly
        // constructed edges
        let verts = Vertex::from_edges(&diff_edges);
        let inter_verts = Vertex::from_edges(&other.outer_edges());

        // The `new_edges` already form the intersection convex polygon
        for vert in verts.iter() {
            for other_vert in inter_verts.iter() {
                if vert == other_vert {
                    continue;
                }
                if colinear_blacklist
                    .get(&vert.idx)
                    .map_or(false, |set| set.contains(&other_vert.idx))
                {
                    continue;
                }
                if split_edges.contains(&Edge::new(*vert, *other_vert)) {
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
                    split_edges.push(Edge::new(*vert, *other_vert));
                    diff_edges.push(Edge::new(*vert, *other_vert));
                }
            }
        }

        trace!(
            "Diff edges after connecting non-intersecting: {:?}",
            diff_edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Result<Vec<_>>>()?
        );

        // Include boundary edges from intersection edges
        let diff_verts = Vertex::from_edges(&diff_edges);
        for edge in other.outer_edges().iter() {
            if diff_verts.contains(&edge.0) && diff_verts.contains(&edge.1) {
                trace!("Trying to include edge {:?} from intersection perimeter",
                    IdxEdge::try_from(edge.clone()));

                // Skip edge that is on the triangle perimeter
                let mut skip = false;
                for tri_edge in self.edges().iter() {
                    let idx_tri_edge = IdxEdge::try_from(*tri_edge).context("Extracting indices for triangle edge")?;
                    let v0_matched = if let Some(IdxIntersection(from_edge_a, from_edge_b)) = edge.0.from {
                        from_edge_a == idx_tri_edge || from_edge_b == idx_tri_edge || edge.0.idx == tri_edge.0.idx || edge.0.idx==tri_edge.1.idx

                    } else {
                        edge.0.idx == tri_edge.0.idx || edge.0.idx==tri_edge.1.idx
                    };
                    let v1_matched = if let Some(IdxIntersection(from_edge_a, from_edge_b)) = edge.1.from {
                        from_edge_a == idx_tri_edge || from_edge_b == idx_tri_edge || edge.1.idx == tri_edge.0.idx || edge.1.idx==tri_edge.1.idx

                    } else {
                        edge.1.idx == tri_edge.0.idx || edge.1.idx==tri_edge.1.idx
                    };
                    if v0_matched && v1_matched {
                        skip = true;
                        trace!("Skip edge {:?} is part of the triangle perimeter and we avoid including the intersection",
                            IdxEdge::try_from(edge.clone()));
                    }
                }
                if skip {
                    continue;
                }

                trace!("Include edge {:?} from intersection perimeter", IdxEdge::try_from(edge.clone()));

                if !diff_edges.contains(edge) {
                    diff_edges.push(*edge);
                } else {
                    warn!(
                        "Had not expect to find an edge from intersection already included in the difference set"
                    );
                }
            }
        }

        trace!(
            "Diff edges to form triangles: {:?}",
            diff_edges
                .iter()
                .map(|e| IdxEdge::try_from(e.clone()))
                .collect::<Vec<_>>()
        );

        let mut tris_diff_vec = self
            .from_edges(diff_edges.clone())
            .context(format!("Construct triangles from the diff edges {:?}", diff_edges))?;

        // In the off chance of the intersection to be a triangle inside self, exclude the triangle
        // formed that matches the intersection
        if other.outer.len() == 3 && other.edges.len() == 3 {
            let polygon = Polygon::try_from(other.outer_edges()).context("Constructing polygon from intersection perimeter")?;
            let other_tri = IdxTriangle::try_from(polygon).unwrap();
            if let Some(pos) = tris_diff_vec.iter().position(|tri| *tri == other_tri) {
                trace!("Excluding triangle {:?} from split result since it overlaps with the intersection triangle {:?}",
                    tris_diff_vec[pos], other_tri);
                tris_diff_vec.remove(pos);
            }
        }

        info!(
            "Diff triangles: {:?}",
            tris_diff_vec
                .iter()
                .map(|tri| tri.verts.map(|v| *v.idx))
                .collect::<Vec<_>>()
        );

        Ok((tris_diff_vec, other))
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
    fn intersect(&self, edges: &Vec<Edge>) -> Result<Vec<Edge>> {
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
            .filter(|Edge(v0, v1)| tri.overlap(&v0.value) && tri.overlap(&v1.value))
            .cloned()
            .collect();
        trace!("Segments after intersection with triangle: {:?}", edges);

        Ok(edges)
    }

    fn split(&self, _edges: Vec<Edge>) -> Result<(Vec<Edge>, Vec<Edge>)> {
        panic!("Not implemented yet");
    }
}

#[cfg(test)]
mod tests {
    // We implement the transformable for the triangle primitive, so we shall use this for tests.
    use super::*;
    use std::sync::Once;

    static INIT: Once = Once::new();

    fn init_logger() {
        INIT.call_once(|| {
            env_logger::Builder::from_default_env()
                .is_test(true) // don't add timestamps
                .try_init()
                .ok();
        });
    }

    #[test]
    fn test_overlaping_triangles() {
        let tri1 = IdxTriangle::new(
            vec![
                Vertex { value: Point3::new(0., 0., 0.), idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(1., 0., 0.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(1., 1., 0.), idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0));
        let tri2 = IdxTriangle::new(
            vec![
                Vertex { value: Point3::new(0., 1., 0.), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(1., 0., 0.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(1., 1., 0.), idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(1));

        assert!(tri1.overlap(&tri2));
    }

    #[test]
    fn test_non_overlaping_triangles() {
        let tri1 = IdxTriangle::new(
            vec![
                Vertex { value: Point3::new(0., 0., 0.), idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(1., 0., 0.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(1., 1., 0.), idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0));

        // Not coplanar
        let tri2 = IdxTriangle::new(
            vec![
                Vertex { value: Point3::new(0., 1., 1.), idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(1., 0., 1.), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(1., 1., 1.), idx: PrimitiveIdx::Global(5), from: None, },
            ], None, PrimitiveIdx::Global(1));
        assert!(!tri1.overlap(&tri2));

        // Coplanar but not overlapping
        let tri2 = IdxTriangle::new(
            vec![
                Vertex { value: Point3::new(0., 0.1, 0.), idx: PrimitiveIdx::Global(6), from: None, },
                Vertex { value: Point3::new(0., 1., 0.), idx: PrimitiveIdx::Global(7), from: None, },
                Vertex { value: Point3::new(0.9, 1., 0.), idx: PrimitiveIdx::Global(8), from: None, },
            ], None, PrimitiveIdx::Global(2));
        assert!(!tri1.overlap(&tri2));

        // Intersecting, but not overlapping in the same plane
        let tri2 = IdxTriangle::new(
            vec![
                Vertex { value: Point3::new(0., 0., -1.), idx: PrimitiveIdx::Global(9), from: None, },
                Vertex { value: Point3::new(1., 0., 1.), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(1., 1., 1.), idx: PrimitiveIdx::Global(5), from: None, },
            ], None, PrimitiveIdx::Global(3));
        assert!(!tri1.overlap(&tri2));
    }

    #[test]
    fn test_overlaping_triangles_one_vertex() -> Result<()> {
        init_logger();
        let tri1 = IdxTriangle::new(
            vec![
                Vertex { value: Point3::new(0., 0., 0.), idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(1., 0., 0.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(1., 1., 0.), idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0),
        );
        let tri2 = IdxTriangle::new(
            vec![
                Vertex { value: Point3::new(0.8, 0.5, 0.), idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(2.,  0.,  0.), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(2.,  1.,  0.), idx: PrimitiveIdx::Global(5), from: None, },
            ], None, PrimitiveIdx::Global(1),
        );

        assert!(tri1.overlap(&tri2));

        info!("Intersecting {:?} and {:?}", tri1.idx, tri2.idx);
        let inter = tri1.intersect(&tri2)?;
        println!("Intersections: {:?}", inter);
        let (new_tri1, inter) = tri1.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();

        println!("New Triangles: {}", new_tri1.len());
        assert_eq!(new_tri1.len(), 4);
        assert_eq!(inter_tris.len(), 1);
        Ok(())
    }

    // Test all triangles split choices
    // 1. resolved => (1,1)
    // 2. edge colinear:
    //   a) 2 edges colinear
    //     i) One edge is the same => (2, 1)
    //     ii) Edges are colinear, but not the same => (3,1)
    //   b) vertex out => (1,1)
    //   c) vertex in:
    //     i) v1=u1 and v2=u2 => (3,1)
    //     ii) v1 = u1, v2 in (u1, u2) => (4,1)
    //     iii) {v1,v2} in (u1, u2) => (5,1)
    // 3. vertex in = 1 (v1)
    //   a) [v1,v2] and [v1,v3] intersects same U edge => (5,3)
    //   b) [v1,v2] and [v1, v2] intersects distinct U edges => (5, 5)
    //   c) One edge of V tri contains U vertex => (4, 3)
    //   d) Two edges of V tri contains U vertices => (3, 3)
    // 4. vertex in = 2 (v1,v2)
    //   a) [v1,v3] and [v2,v3] intersects U edges
    //     i) same U edge => (7,3)
    //     ii) 2 distinct U edges => (7, 5)
    //   b) One edge of V contains U vertex => (6, 3)
    // 5. vertex in = 3 => (7,1)
    // 6. vertex in = 0
    //   a) 2 edges of V each intersect with any 2 edges of U
    //     i) Vertex of U inside V => (5, 7)
    //     ii) Vertex of U on edge of V => (5, 6)
    //   b) 1 edge of V intersects with 2 edges of U
    //     i) One vertex of U on edge of V => (3, 6)
    //     ii) Two vertices of U on edges of V => (3, 5)
    //     iii) Two vertices of U inside V => (3, 7)
    //   c) 3 edges of V intersect 3 edges of U => (7,7)

    fn count_colinear(us: &[Segment; 3], vs: &[Segment; 3]) -> usize {
        let mut cnt = 0;
        for u in us {
            for v in vs {
                if u.colinear(&v) {
                    cnt += 1;
                }
            }
        }
        assert!(
            cnt <= 3,
            "Can't have more than 3 colinear cross edges between 2 triangles, but counted {}",
            cnt
        );
        cnt
    }

    fn count_verts_in(u_tri: &IdxTriangle, v_tri: &IdxTriangle) -> usize {
        v_tri.verts.iter().fold(0, |cnt, v| {
            if Into::<Triangle>::into(u_tri.clone()).vertex_in(v.value) {
                cnt + 1
            } else {
                cnt
            }
        })
    }

    #[test]
    fn test_matching_triangles_split() -> Result<()> {
        let tri1 = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0., 0.), idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 1., 0.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(0., 1., 1.), idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0));
        let mut tri2 = tri1.clone();
        tri2.idx = PrimitiveIdx::Global(1);

        assert!(tri1.overlap(&tri2));
        let tri1_segs: [Segment; 3] = tri1.edges().map(|e| e.into());
        let tri2_segs: [Segment; 3] = tri2.edges().map(|e| e.into());
        assert_eq!(count_colinear(&tri1_segs, &tri2_segs), 3);

        let inter = tri1.intersect(&tri2)?;
        let (new_tri1, inter) = tri1.split(inter)?;
        let (new_tri2, inter) = tri2.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();

        assert_eq!(new_tri1.len(), 0);
        assert_eq!(new_tri2.len(), 0);
        assert_eq!(inter_tris.len(), 1);

        assert_eq!(inter_tris[0].tri(), tri1.tri());
        Ok(())
    }

    #[test]
    fn test_colinear_edges() -> Result<()> {
        init_logger();
        let tri_u = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0., 1.),  idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 0., -1.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(0., 1., 0.),  idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0));
        let segs_u = tri_u.segs();

        // a) i)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0.,  1.),  idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 0.,  -1.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(0., 0.5, 0.5), idx: PrimitiveIdx::Global(3), from: None, },
            ], None, PrimitiveIdx::Global(1));
        let segs_v = tri_v.segs();
        assert_eq!(count_colinear(&segs_u, &segs_v), 2);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 1);
        assert_eq!(new_tris_v.len(), 0);
        assert_eq!(inter_tris.len(), 1);

        // a) ii)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0.,  1.),  idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 0.,  0.),  idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(0., 0.5, 0.5), idx: PrimitiveIdx::Global(3), from: None, },
            ], None, PrimitiveIdx::Global(2));
        let segs_v = tri_v.segs();
        assert_eq!(count_colinear(&segs_u, &segs_v), 2);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 2);
        assert_eq!(new_tris_v.len(), 0);
        assert_eq!(inter_tris.len(), 1);

        // b)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0.,  1.),   idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 0.,  -0.8), idx: PrimitiveIdx::Global(5), from: None, },
                Vertex { value: Point3::new(0., -1., 0.),   idx: PrimitiveIdx::Global(6), from: None, },
            ], None, PrimitiveIdx::Global(3));
        let segs_v = tri_v.segs();
        assert_eq!(count_colinear(&segs_u, &segs_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 1);
        assert_eq!(new_tris_v.len(), 1);
        assert_eq!(inter_tris.len(), 0);

        // c) vert in
        // i) // WARN: Duplicate of vert_in=1 d)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0.,  1.),  idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 0.,  -1.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(0., 0.5, 0.),  idx: PrimitiveIdx::Global(8), from: None, },
            ], None, PrimitiveIdx::Global(4));
        let segs_v = tri_v.segs();
        assert_eq!(count_colinear(&segs_u, &segs_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 2);
        assert_eq!(new_tris_v.len(), 0);
        assert_eq!(inter_tris.len(), 1);

        // ii) // WARN: Duplicate of vert_in=1 c)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0.,  1.),   idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 0.,  -0.5), idx: PrimitiveIdx::Global(9), from: None, },
                Vertex { value: Point3::new(0., 0.5, 0.),   idx: PrimitiveIdx::Global(8), from: None, },
            ], None, PrimitiveIdx::Global(5));
        let segs_v = tri_v.segs();
        assert_eq!(count_colinear(&segs_u, &segs_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        println!("Intersection: {:?}", inter.edges.iter().map(|e| IdxEdge::try_from(e.clone())).collect::<Vec<_>>());
        println!("New tris U: {:?}", new_tris_u.iter().map(|tri| tri.verts.map(|v| v.idx)).collect::<Vec<_>>());
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 3);
        assert_eq!(new_tris_v.len(), 0);
        assert_eq!(inter_tris.len(), 1);

        // iii) // WARN: Duplicate of vert_in=1 b)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0.,  0.5),  idx: PrimitiveIdx::Global(10), from: None, },
                Vertex { value: Point3::new(0., 0.,  -0.5), idx: PrimitiveIdx::Global(9),  from: None, },
                Vertex { value: Point3::new(0., 0.5, 0.),   idx: PrimitiveIdx::Global(8),  from: None, },
            ], None, PrimitiveIdx::Global(6));
        let segs_v = tri_v.segs();
        assert_eq!(count_colinear(&segs_u, &segs_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 4);
        assert_eq!(new_tris_v.len(), 0);
        assert_eq!(inter_tris.len(), 1);

        // iv) All  vertices of V on edges of U
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0.,  0.5),  idx: PrimitiveIdx::Global(10), from: None, },
                Vertex { value: Point3::new(0., 0.,  -0.5), idx: PrimitiveIdx::Global(9),  from: None, },
                Vertex { value: Point3::new(0., 0.5, 0.5),  idx: PrimitiveIdx::Global(11), from: None, },
            ], None, PrimitiveIdx::Global(7));
        let segs_v = tri_v.segs();
        assert_eq!(count_colinear(&segs_u, &segs_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 3);
        assert_eq!(new_tris_v.len(), 0);
        assert_eq!(inter_tris.len(), 1);
        Ok(())
    }

    #[test]
    fn test_one_vert_in() -> Result<()> {
        let tri_u = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 0.5, 0.5),  idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 0.5, -0.5), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(0., 1.5, 0.),   idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0));

        // a)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 1.,  0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0.,  0.5), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(0., 0., -0.5), idx: PrimitiveIdx::Global(5), from: None, },
            ], None, PrimitiveIdx::Global(1));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 4);
        assert_eq!(new_tris_v.len(), 2);
        assert_eq!(inter_tris.len(), 1);

        // b)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 1.,  0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 2.,  0.5), idx: PrimitiveIdx::Global(6), from: None, },
                Vertex { value: Point3::new(0., 2., -0.5), idx: PrimitiveIdx::Global(7), from: None, },
            ], None, PrimitiveIdx::Global(2));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 3);
        assert_eq!(new_tris_v.len(), 3);
        assert_eq!(inter_tris.len(), 2);

        // c)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 1.,  0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0.,  0.5), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(0., 0., -1.), idx: PrimitiveIdx::Global(8), from: None, },
            ], None, PrimitiveIdx::Global(3));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 3);
        assert_eq!(new_tris_v.len(), 2);
        assert_eq!(inter_tris.len(), 1);

        // d)
        let tri_v = IdxTriangle::new( vec![
                Vertex { value: Point3::new(0., 1.,  0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0.,  1.), idx: PrimitiveIdx::Global(9), from: None, },
                Vertex { value: Point3::new(0., 0., -1.), idx: PrimitiveIdx::Global(8), from: None, },
            ], None, PrimitiveIdx::Global(4));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 1);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 2);
        assert_eq!(new_tris_v.len(), 2);
        assert_eq!(inter_tris.len(), 1);
        Ok(())
    }

    #[test]
    fn test_two_verts_in() -> Result<()> {
        let tri_u = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 2., 0.),  idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 0., -1.), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(0., 0., 1.),   idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0));

        // a) i)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 0.5, 0.5),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0.5, -0.5), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(0., -1., 0.),   idx: PrimitiveIdx::Global(5), from: None, },
            ], None, PrimitiveIdx::Global(1));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 2);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 5);
        assert_eq!(new_tris_v.len(), 1);
        assert_eq!(inter_tris.len(), 2);

        // a) ii)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 0.5, 0.5),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0.5, -0.5), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(0., 3., 0.),   idx: PrimitiveIdx::Global(6), from: None, },
            ], None, PrimitiveIdx::Global(2));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 2);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 4);
        assert_eq!(new_tris_v.len(), 2);
        assert_eq!(inter_tris.len(), 3);

        // b)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 0.5, 0.5),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0.5, 0.), idx: PrimitiveIdx::Global(7), from: None, },
                Vertex { value: Point3::new(0., 3., 0.),   idx: PrimitiveIdx::Global(6), from: None, },
            ], None, PrimitiveIdx::Global(3));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 2);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 4);
        assert_eq!(new_tris_v.len(), 1);
        assert_eq!(inter_tris.len(), 2);
        Ok(())
    }

    #[test]
    fn test_all_verts_in() -> Result<()> {
        init_logger();
        let tri_u = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., -1., 0.5),  idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., -1., -0.5), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(0., 1., 0.),   idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0));
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., -0.5, 0.25),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., -0.5, -0.25), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(0., 0.5, 0.),   idx: PrimitiveIdx::Global(5), from: None, },
            ], None, PrimitiveIdx::Global(1));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 3);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 6);
        assert_eq!(new_tris_v.len(), 0);
        assert_eq!(inter_tris.len(), 1);
        Ok(())
    }

    #[test]
    fn test_verts_out() -> Result<()> {
        init_logger();
        let tri_u = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 1., 0.5),  idx: PrimitiveIdx::Global(0), from: None, },
                Vertex { value: Point3::new(0., 1., -0.5), idx: PrimitiveIdx::Global(1), from: None, },
                Vertex { value: Point3::new(0., -1., 0.),   idx: PrimitiveIdx::Global(2), from: None, },
            ], None, PrimitiveIdx::Global(0));

        // a) i)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 2., 0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0., -0.4), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(0., 0., 1.5),   idx: PrimitiveIdx::Global(5), from: None, },
            ], None, PrimitiveIdx::Global(1));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 0);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 2);
        assert_eq!(new_tris_v.len(), 4);
        assert_eq!(inter_tris.len(), 3);

        // a) ii)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 2., 0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0., -0.4), idx: PrimitiveIdx::Global(4), from: None, },
                Vertex { value: Point3::new(0., 0., 1.),   idx: PrimitiveIdx::Global(6), from: None, },
            ], None, PrimitiveIdx::Global(2));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 0);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 2);
        assert_eq!(new_tris_v.len(), 3);
        assert_eq!(inter_tris.len(), 3);

        // b) i)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 2., 0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0., -1.), idx: PrimitiveIdx::Global(7), from: None, },
                Vertex { value: Point3::new(0., 0., 1.5),   idx: PrimitiveIdx::Global(8), from: None, },
            ], None, PrimitiveIdx::Global(3));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 0);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 1);
        assert_eq!(new_tris_v.len(), 4);
        assert_eq!(inter_tris.len(), 2);

        // b) ii)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 2., 0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0., -1.), idx: PrimitiveIdx::Global(7), from: None, },
                Vertex { value: Point3::new(0., 0., 1.),   idx: PrimitiveIdx::Global(6), from: None, },
            ], None, PrimitiveIdx::Global(4));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 0);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 1);
        assert_eq!(new_tris_v.len(), 3);
        assert_eq!(inter_tris.len(), 2);

        // b) iii)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 2., 0.),  idx: PrimitiveIdx::Global(3), from: None, },
                Vertex { value: Point3::new(0., 0., -1.5), idx: PrimitiveIdx::Global(9), from: None, },
                Vertex { value: Point3::new(0., 0., 1.5),   idx: PrimitiveIdx::Global(8), from: None, },
            ], None, PrimitiveIdx::Global(5));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 0);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 1);
        assert_eq!(new_tris_v.len(), 5);
        assert_eq!(inter_tris.len(), 2);

        // c)
        let tri_v = IdxTriangle::new(vec![
                Vertex { value: Point3::new(0., 1.1, 0.),  idx: PrimitiveIdx::Global(10), from: None, },
                Vertex { value: Point3::new(0., 0., -1.), idx: PrimitiveIdx::Global(7), from: None, },
                Vertex { value: Point3::new(0., 0., 1.),   idx: PrimitiveIdx::Global(6), from: None, },
            ], None, PrimitiveIdx::Global(6));
        assert_eq!(count_verts_in(&tri_u, &tri_v), 0);
        let inter = tri_u.intersect(&tri_v)?;
        let (new_tris_u, inter) = tri_u.split(inter)?;
        let (new_tris_v, inter) = tri_v.split(inter)?;
        let inter_tris: Vec<IdxTriangle> = inter.try_into().unwrap();
        assert_eq!(new_tris_u.len(), 3);
        assert_eq!(new_tris_v.len(), 3);
        assert_eq!(inter_tris.len(), 4);
        Ok(())
    }

    //#[test]
    //fn test_segments_split() {
    //    let tri_u = Triangle::new([
    //        Point3::new(0., 2., 0.),
    //        Point3::new(0., 0., -1.),
    //        Point3::new(0., 0., 1.),
    //    ]);

    //    // a) ii)
    //    let tri_v = Triangle::new([
    //        Point3::new(0., 0.5, 0.5),
    //        Point3::new(0., 0.5, -0.5),
    //        Point3::new(0., 3., 0.),
    //    ]);
    //    assert_eq!(count_verts_in(&tri_u, &tri_v), 2);

    //    let (tris_u, split_segments) = tri_u.split_transparent(&tri_v);
    //    assert_eq!(tris_u.len(), 7);

    //    println!("Split segments: {:?}", split_segments);

    //    let tris_v = tri_v.split(&split_segments);
    //    assert_eq!(tris_v.len(), 5);
    //}
}
