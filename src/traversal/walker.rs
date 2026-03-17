//! Module for walking the [`Smiles`] graph and tracking progress

use std::collections::HashSet;

use crate::{
    errors::SmilesError, smiles::Smiles,
    traversal::visitor_trait::Visitor,
};

/// Traverses the [`Smiles`] graph via a depth first search
pub fn walk<V: Visitor>(smiles: &Smiles, visitor: &mut V) -> Result<(), SmilesError> {
    let mut visited_nodes: HashSet<usize> = HashSet::new();

    for node in smiles.nodes() {
        if !visited_nodes.contains(&node.id()) {
            dfs(smiles, visitor, node.id(), &mut visited_nodes)?;
        }
    }

    Ok(())
}
fn dfs<V: Visitor>(
    smiles: &Smiles,
    visitor: &mut V,
    current_id: usize,
    visited_nodes: &mut HashSet<usize>,
) -> Result<(), SmilesError> {
    if visited_nodes.contains(&current_id) {
        return Ok(());
    }
    visited_nodes.insert(current_id);
    visitor.enter_node(smiles, current_id)?;
    let bonds = smiles.edges_for_node(current_id);
    for bond in bonds {
        if let Some(other_id) = bond.other(current_id) {
            if !visited_nodes.contains(&other_id) {
                visitor.tree_edge(smiles, *bond)?;
                dfs(smiles, visitor, other_id, visited_nodes)?;
            }
        }
    }
    visitor.exit_node(smiles, current_id)?;
    Ok(())
}
