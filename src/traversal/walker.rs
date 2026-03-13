//! Module for walking the [`Smiles`] graph and tracking progress

use std::collections::HashSet;

use crate::{errors::SmilesError, smiles::Smiles, traversal::visitor_trait::Visitor};

/// Traverses the [`Smiles`] graph via a depth first search
pub fn walk<V: Visitor>(smiles: &Smiles, visitor: &mut V) -> Result<(), SmilesError> {
    let mut visited_nodes: HashSet<usize> = HashSet::new();

    for node in smiles.nodes() {
        visitor.enter_node(smiles, node.id())?;
        visited_nodes.insert(node.id());

        

        visitor.exit_node(smiles, node.id())?;
    }
    Ok(())
}
