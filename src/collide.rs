//! Collide trait.
/// As this trait is geometry dependent, we will test the collide implementation per implementing type.
pub trait Collide<T> {
    /// Check for an overlapping collision.
    fn overlap(&self, other: &T) -> bool;
}

