//! Module for graph traversal visitors.

use crate::{bond::Bond, errors::SmilesError, smiles::Smiles};

/// Trait that defines visitors to [`Smiles`] nodes
pub trait Visitor {
    fn enter_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError>;
    fn exit_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError>;

    fn tree_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError>;

    fn cycle_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError>;
}
