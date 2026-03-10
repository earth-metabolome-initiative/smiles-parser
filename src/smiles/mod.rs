//! Represents a SMILES structure.


use crate::{
    atom::atom_node::AtomNode,
    bond::{Bond, bond_edge::BondEdge},
    errors::SmilesError,
};

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
    ///
    /// # Errors
    /// - Returns a [`SmilesError::NodeIdInvalid`] if a node cannot be found in
    ///   the edge list
    pub fn push_edge(
        &mut self,
        node_a: usize,
        node_b: usize,
        bond: Bond,
    ) -> Result<(), SmilesError> {
        self.atom_nodes.sort();
        // use the NodeIdInvalid for err
        if self
            .atom_nodes
            .binary_search_by_key(&node_a, super::atom::atom_node::AtomNode::id)
            .is_err()
        {
            return Err(SmilesError::NodeIdInvalid(node_a));
        }
        if self
            .atom_nodes
            .binary_search_by_key(&node_b, super::atom::atom_node::AtomNode::id)
            .is_err()
        {
            return Err(SmilesError::NodeIdInvalid(node_b));
        }
        let bond_edge = BondEdge::new(node_a, node_b, bond);
        self.bond_edges.push(bond_edge);

        Ok(())
    }
    /// Returns slice of the nodes
    #[must_use]
    pub fn nodes(&self) -> &[AtomNode] {
        &self.atom_nodes
    }
    /// Returns mutable slice of nodes
    #[must_use]
    pub fn nodes_mut(&mut self) -> &mut [AtomNode] {
        &mut self.atom_nodes
    }
    /// Returns a node if exists from an `id`
    #[must_use]
    pub fn node_by_id(&self, id: usize) -> Option<&AtomNode> {
        self.atom_nodes.iter().find(|node| node.id() == id)
    }
    /// Returns slice of the edges
    #[must_use]
    pub fn edges(&self) -> &[BondEdge] {
        &self.bond_edges
    }
    /// Returns mutable slice of edges
    pub fn edges_mut(&mut self) -> &mut [BondEdge] {
        &mut self.bond_edges
    }
    /// Find all neighboring edges bonded to the current node id
    #[must_use]
    pub fn neighbors(&self, id: usize) -> Vec<(usize, Bond)> {
        self.bond_edges
            .iter()
            .filter_map(|edge| {
                if edge.node_a() == id {
                    Some((edge.node_b(), edge.bond().to_owned()))
                } else if edge.node_b() == id {
                    Some((edge.node_a(), edge.bond().to_owned()))
                } else {
                    None
                }
            })
            .collect()
    }
}

impl Default for Smiles {
    fn default() -> Self {
        Self::new()
    }
}

