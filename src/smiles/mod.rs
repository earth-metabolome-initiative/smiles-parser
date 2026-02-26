//! Represents a SMILES structure.

use crate::{atom::atom_node::AtomNode, bond::Bond, bond::bond_edge::BondEdge};

mod from_str;

/// Represents a SMILES structure.
pub struct Smiles {
    atom_nodes: Vec<AtomNode>,
    bond_edges: Vec<BondEdge>,
}

impl Smiles {
    /// creates a new instance
    #[must_use]
    pub fn new() -> Self {
        Self { atom_nodes: Vec::new(), bond_edges: Vec::new() }
    }
    /// Pushes an AtomNode
    pub fn push_node(&mut self, node: AtomNode) {
        self.atom_nodes.push(node);
    }
    /// adds an edge from two nodes and the [`Bond`]
    pub fn push_edge(&mut self, node_a: usize, node_b: usize, bond: Bond) {
        let bond_edge = BondEdge::new(node_a, node_b, bond);
        self.bond_edges.push(bond_edge);
    }
}

impl Default for Smiles {
    fn default() -> Self { Self::new() }
}