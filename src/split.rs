use crate::collide::Collide;

pub trait Split<T, P>: Collide<T>
where
    Self: Sized
{
    type Inst;
    // Using the intersection primitives, split object that forms the difference between the
    // original object and the one intersected with
    fn split(&self, other: &P) -> (Self::Inst, P);

    // Find the intersection between two objects and return the primitives that described the
    // intersection
    fn intersect(&self, other: &T) -> P;
}
