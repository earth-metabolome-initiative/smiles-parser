//! Module for bonds as edges for a graph structure

use crate::bond::Bond;

/// Contains the two ID's of the [`AtomNode`] that are connected via the
/// [`Bond`]
pub struct BondEdge {
    /// The first node
    node_a: usize,
    /// The Second node
    node_b: usize,
    /// The bond between the nodes
    bond: Bond,
}

impl BondEdge {
    /// Creates a new edge
    #[must_use]
    pub fn new(node_a: usize, node_b: usize, bond: Bond) -> Self {
        Self { node_a, node_b, bond }
    }
    /// Returns the specified [`Bond`]
    #[must_use]
    pub fn bond(&self) -> &Bond {
        &self.bond
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
}
