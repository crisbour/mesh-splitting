mod segment;
mod collide;
mod split;
mod save;
mod triangle;
mod idx_triangle;
pub mod primitives;
pub mod mesh;
pub mod export;

use std::ops::Add;

use nalgebra::{Point3, Vector3};
pub use collide::Collide;
pub use save::Save;
pub use split::Split;
pub use segment::{Segment, SplitEdges};
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
    pub fn null() -> Self {
        Self {
            mins: Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            maxs: Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }
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

impl Add<&Point3<f64>> for Aabb {
    type Output = Self;
    fn add(self, rhs: &Point3<f64>) -> Self::Output {
        self + *rhs
    }
}

impl Add<Point3<f64>> for Aabb {
    type Output = Self;
    fn add(self, rhs: Point3<f64>) -> Self::Output {
        //let mins = self.mins.simd_min(rhs);
        let mins = Point3::new(self.mins[0].min(rhs[0]), self.mins[1].min(rhs[1]), self.mins[2].min(rhs[2]));
        //let maxs = self.maxs.simd_max(rhs);
        let maxs = Point3::new(self.maxs[0].max(rhs[0]), self.maxs[1].max(rhs[1]), self.maxs[2].max(rhs[2]));
        Self { mins, maxs }
    }
}

impl Add<Aabb> for Aabb {
    type Output = Self;
    fn add(self, rhs: Aabb) -> Self::Output {
        self + rhs.mins + rhs.maxs
    }
}

impl Collide<Aabb> for Aabb {
    fn overlap(&self, aabb: &Aabb) -> bool {
        self.mins <= aabb.maxs && self.maxs >= aabb.mins
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aabb_union() {
        let aabb1 = Aabb {
            mins: Point3::new(0.0, 0.0, 0.0),
            maxs: Point3::new(1.0, 1.0, 1.0),
        };
        let aabb2 = Aabb {
            mins: Point3::new(2.0, 2.0, 2.0),
            maxs: Point3::new(3.0, 3.0, 3.0),
        };
        let union = aabb1 + aabb2;
        assert_eq!(union.mins, Point3::new(0.0, 0.0, 0.0));
        assert_eq!(union.maxs, Point3::new(3.0, 3.0, 3.0));
    }
}
