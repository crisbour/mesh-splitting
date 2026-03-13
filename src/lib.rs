mod segment;
mod collide;
mod split;
mod triangle;
mod idx_triangle;
pub mod primitives;
pub mod mesh;

use nalgebra::{Point3, Vector3};
pub use segment::{Segment, SplitEdges};
pub use collide::Collide;
pub use split::Split;
pub use triangle::Triangle;
pub use idx_triangle::IdxTriangle;


#[derive(Clone, Debug)]
pub struct Aabb {
    /// Minimum bound.
    mins: Point3<f64>,
    /// Maximum bound.
    maxs: Point3<f64>,
}

impl Aabb {
    /// Calculate the widths.
    pub fn widths(&self) -> Vector3<f64> {
        self.maxs - self.mins
    }

    /// Calculate the half-widths.
    pub fn half_widths(&self) -> Vector3<f64> {
        self.widths() * 0.5
    }

    pub fn centre(&self) -> Point3<f64> {
        nalgebra::center(&self.mins, &self.maxs).into()
    }
}

impl Collide<Aabb> for Aabb {
    fn overlap(&self, aabb: &Aabb) -> bool {
        self.mins <= aabb.maxs && self.maxs >= aabb.mins
    }
}

