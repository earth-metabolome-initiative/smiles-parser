//! Module for graph traversal visitors.

use crate::{
    bond::{Bond, bond_edge::BondEdge},
    errors::SmilesError,
    smiles::Smiles,
};

/// Trait that defines visitors to [`Smiles`] nodes
pub trait Visitor {
    /// Called when traversal first enters a node
    ///
    /// # Errors
    /// - Will return a [`SmilesError`] with context of error if fails to enter
    ///   a node
    fn enter_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError>;
    /// Called when a traversal finishes processing a node.
    ///
    /// # Errors
    /// - Will return a [`SmilesError`] with context of error if fails to exit a
    ///   node
    fn exit_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError>;
    /// Called when traversal follows an edge to an unvisited
    /// node
    ///
    /// # Errors
    /// - Will return a [`SmilesError`] with context of error if fails to
    ///   evaluate a possible tree edge
    fn tree_edge(&mut self, smiles: &Smiles, bond_edge: BondEdge) -> Result<(), SmilesError>;
    /// Called when traversal encounters an edge that closes a
    /// cycle (a ring)
    ///
    /// # Errors
    /// - Will return a [`SmilesError`] with context of error if fails to
    ///   evaluate a possible cycle edge
    fn cycle_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError>;
    /// Called before traversing a side branch has begun
    ///
    /// # Errors
    /// - Will return a [`SmilesError`] with context of error if fails to open a
    ///   branch
    fn open_branch(&mut self, smiles: &Smiles, from: usize, to: usize) -> Result<(), SmilesError>;
    /// Called after traversing a side branch has concluded
    ///
    /// # Errors
    /// - Will return a [`SmilesError`] with context of error if fails to close
    ///   the branch
    fn close_branch(&mut self, smiles: &Smiles, from: usize, to: usize) -> Result<(), SmilesError>;
    /// Called before entering a new sub graph
    ///
    /// # Errors
    /// - Will return a [`SmilesError`] with context of error if fails to start
    ///   a component
    fn start_component(
        &mut self,
        smiles: &Smiles,
        root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError>;
    /// Called after exiting a sub graph
    ///
    /// # Errors
    /// - Will return a [`SmilesError`] with context of error if fails to finish
    ///   a component
    fn finish_component(
        &mut self,
        smiles: &Smiles,
        root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError>;
}
