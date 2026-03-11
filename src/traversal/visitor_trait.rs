//! Module for graph traversal visitors.

use crate::{bond::Bond, errors::SmilesError, smiles::Smiles};

/// Trait that defines visitors to [`Smiles`] nodes
pub trait Visitor {
    /// Intended to be called when traversal first enters a node
    fn enter_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError>;
    /// Intended to be called when a traversal finishes processing a node.
    fn exit_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError>;
    /// Intended to be called when traversal follows an edge to an unvisited node
    fn tree_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError>;
    /// Intended to be called when traversal encounters an edge that closes a cycle (a ring)
    fn cycle_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError>;
}
