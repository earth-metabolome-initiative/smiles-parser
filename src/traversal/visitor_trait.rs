//! Module for graph traversal visitors.

use crate::{bond::Bond, errors::SmilesError, smiles::Smiles};

/// Trait that defines visitors to [`Smiles`] nodes
pub trait Visitor {
    /// Called when traversal first enters a node
    fn enter_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError>;
    /// Called when a traversal finishes processing a node.
    fn exit_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError>;
    /// Called when traversal follows an edge to an unvisited
    /// node
    fn tree_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError>;
    /// Called when traversal encounters an edge that closes a
    /// cycle (a ring)
    fn cycle_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError>;
    /// Called before traversing a side branch has begun
    fn open_branch(&mut self, smiles: &Smiles, from: usize, to: usize) -> Result<(), SmilesError>;
    /// Called after traversing a side branch has concluded
    fn close_branch(&mut self, smiles: &Smiles, from: usize, to: usize) -> Result<(), SmilesError>;
    /// Called before entering a new sub graph
    fn start_component(
        &mut self,
        smiles: &Smiles,
        root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError>;
    /// Called after exiting a sub graph
    fn finish_component(
        &mut self,
        smiles: &Smiles,
        root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError>;
}
