//! Module for walking the [`Smiles`] graph and tracking progress

use std::collections::HashSet;

use crate::{
    atom::atom_node::AtomNode, errors::SmilesError, smiles::Smiles,
    traversal::visitor_trait::Visitor,
};

/// Traverses the [`Smiles`] graph via a depth first search
pub fn walk<V: Visitor>(smiles: &Smiles, visitor: &mut V) -> Result<(), SmilesError> {
    let mut visited_nodes: HashSet<usize> = HashSet::new();
    visitor.start(smiles)?;
    
    dfs(smiles, visitor, &mut visited_nodes);

    Ok(())
}

fn dfs<V: Visitor>(smiles: &Smiles, visitor: &mut V, visited_nodes: &mut HashSet<usize>) {
    for node in smiles.nodes() {
        if !visited_nodes.contains(&node.id()) {
            visited_nodes.insert(node.id());
        }
    }
}
