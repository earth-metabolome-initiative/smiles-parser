//! Module for bonds as edges for a graph structure

use std::fmt;

use crate::bond::{Bond, ring_num::RingNum};

/// Contains the two ID's of the `AtomNode` that are connected via the
/// [`Bond`]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondEdge {
    /// The first node
    node_a: usize,
    /// The Second node
    node_b: usize,
    /// The bond between the nodes
    bond: Bond,
    /// Possible Ring Number for this edge
    ring_num: Option<RingNum>,
}

impl BondEdge {
    /// Creates a new edge
    #[must_use]
    pub fn new(node_a: usize, node_b: usize, bond: Bond, ring_num: Option<RingNum>) -> Self {
        Self { node_a, node_b, bond, ring_num }
    }
    /// Returns the specified [`Bond`]
    #[must_use]
    pub fn bond(&self) -> Bond {
        self.bond
    }
    /// Returns a tuple of the two vertices
    #[must_use]
    pub fn vertices(&self) -> (usize, usize) {
        (self.node_a, self.node_b)
    }
    /// Returns the first vertex
    #[must_use]
    pub fn node_a(&self) -> usize {
        self.node_a
    }
    /// Returns the second vertex
    #[must_use]
    pub fn node_b(&self) -> usize {
        self.node_b
    }
    /// Checks if the `BondEdge` contains the specified node id as `usize`
    #[must_use]
    pub fn contains(&self, node_id: usize) -> bool {
        self.node_a == node_id || self.node_b == node_id
    }
    /// Returns the other node id (if it exists)
    #[must_use]
    pub fn other(&self, node_id: usize) -> Option<usize> {
        if self.node_a == node_id {
            Some(self.node_b)
        } else if self.node_b == node_id {
            Some(self.node_a)
        } else {
            None
        }
    }
    /// Returns the ring number assigned to this edge, if any.
    #[must_use]
    pub fn ring_num(&self) -> Option<RingNum> {
        self.ring_num
    }
    /// Updates the [`RingNum`] value
    pub fn set_ring_num(&mut self, ring_num: Option<RingNum>) {
        self.ring_num = ring_num
    }
    /// returns the [`RingNum`] value
    #[must_use]
    pub fn ring_num_val(&self) -> Option<u8> {
        self.ring_num.map(|num| num.get())
    }
}

impl fmt::Display for BondEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.bond {
            Bond::Single => Ok(()),
            _ => write!(f, "{}", self.bond),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bond::{Bond, bond_edge::BondEdge, ring_num::RingNum};

    #[test]
    fn test_bond_edge_new_and_accessors() {
        let edge = BondEdge::new(3, 7, Bond::Double, None);

        assert_eq!(edge.node_a(), 3);
        assert_eq!(edge.node_b(), 7);
        assert_eq!(edge.vertices(), (3, 7));
        assert_eq!(edge.bond(), Bond::Double);
        assert_eq!(edge.ring_num(), None);
    }

    #[test]
    fn test_bond_edge_fmt_all_arms() {
        let cases = [
            (BondEdge::new(0, 1, Bond::Single, Some(RingNum::try_new(2).unwrap())), ""),
            (BondEdge::new(0, 1, Bond::Double, Some(RingNum::try_new(2).unwrap())), "="),
            (BondEdge::new(0, 1, Bond::Triple, Some(RingNum::try_new(2).unwrap())), "#"),
            (BondEdge::new(0, 1, Bond::Quadruple, Some(RingNum::try_new(2).unwrap())), "$"),
            (BondEdge::new(0, 1, Bond::Aromatic, Some(RingNum::try_new(2).unwrap())), ":"),
            (BondEdge::new(0, 1, Bond::Up, Some(RingNum::try_new(2).unwrap())), "/"),
            (BondEdge::new(0, 1, Bond::Down, Some(RingNum::try_new(2).unwrap())), "\\"),
        ];

        for (edge, expected) in cases {
            assert_eq!(expected, edge.to_string());
        }
    }

    #[test]
    fn test_bond_edge_vertices_preserve_order() {
        let edge = BondEdge::new(10, 2, Bond::Single, None);
        assert_eq!(edge.vertices(), (10, 2));
    }
}
