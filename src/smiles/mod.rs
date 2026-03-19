//! Represents a SMILES structure.

use std::fmt;

use crate::{
    atom::atom_node::AtomNode,
    bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
    errors::SmilesError,
    traversal::{render_visitor::RenderVisitor, walker::walk},
};

mod from_str;

/// Represents a SMILES structure.
#[derive(Debug)]
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
    ///
    /// # Errors
    /// - Returns [`SmilesError::DuplicateNodeId`] if node id already exists
    pub fn push_node(&mut self, node: AtomNode) -> Result<(), SmilesError> {
        if self.node_by_id(node.id()).is_some() {
            return Err(SmilesError::DuplicateNodeId(node.id()));
        }
        self.atom_nodes.push(node);
        Ok(())
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
        ring_num: Option<RingNum>,
    ) -> Result<(), SmilesError> {
        if !self.contains_node_id(node_a) {
            return Err(SmilesError::NodeIdInvalid(node_a));
        }
        if !self.contains_node_id(node_b) {
            return Err(SmilesError::NodeIdInvalid(node_b));
        }
        self.bond_edges.push(BondEdge::new(node_a, node_b, bond, ring_num));
        Ok(())
    }
    /// Returns `bool` for if the `AtomNode` `id` exists in the set
    fn contains_node_id(&self, id: usize) -> bool {
        self.atom_nodes.iter().any(|node| node.id() == id)
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
    /// Returns a normalized edge key with node IDs in ascending order.
    #[must_use]
    pub fn edge_key(a: usize, b: usize) -> (usize, usize) {
        if a < b { (a, b) } else { (b, a) }
    }
    /// Returns the first [`BondEdge`] connecting the given pair of node IDs.
    #[must_use]
    pub fn edge_for_node_pair(&self, nodes: (usize, usize)) -> Option<&BondEdge> {
        let target = Self::edge_key(nodes.0, nodes.1);
        self.bond_edges.iter().find(|b| {
            let (a, c) = b.vertices();
            Self::edge_key(a, c) == target
        })
    }
    /// Returns all edges for a given node `id`
    #[must_use]
    pub fn edges_for_node(&self, id: usize) -> Vec<&BondEdge> {
        self.bond_edges.iter().filter(|b| b.contains(id)).collect()
    }
    /// Returns mutable slice of edges
    pub fn edges_mut(&mut self) -> &mut [BondEdge] {
        &mut self.bond_edges
    }
    /// Renders the `Smiles` graph into a SMILES string
    ///
    /// # Errors
    /// - Returns a [`SmilesError`] if the graph fails to walk
    pub fn render(&self) -> Result<String, SmilesError> {
        let mut render_visitor = RenderVisitor::new();
        walk(self, &mut render_visitor)?;
        Ok(render_visitor.into_string())
    }
}

impl Default for Smiles {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Smiles {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // need to map the smiles error to the `fmt::Error`
        let rendered_smiles = self.render().map_err(|_| fmt::Error)?;
        write!(f, "{rendered_smiles}")
    }
}
