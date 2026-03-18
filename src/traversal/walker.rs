//! Module for walking the [`Smiles`] graph and tracking progress

use std::collections::HashSet;

use crate::{errors::SmilesError, smiles::Smiles, traversal::visitor_trait::Visitor};

/// Traverses the [`Smiles`] graph via a depth first search
pub fn walk<V: Visitor>(smiles: &Smiles, visitor: &mut V) -> Result<(), SmilesError> {
    let mut visited_nodes: HashSet<usize> = HashSet::new();
    let mut visited_edges: HashSet<(usize, usize)> = HashSet::new();
    let mut component_index = 0;

    for node in smiles.nodes() {
        if !visited_nodes.contains(&node.id()) {
            visitor.start_component(smiles, node.id(), component_index)?;
            dfs(smiles, visitor, node.id(), &mut visited_nodes, &mut visited_edges)?;
            visitor.finish_component(smiles, node.id(), component_index)?;
            component_index += 1;
        }
    }

    Ok(())
}

fn dfs<V: Visitor>(
    smiles: &Smiles,
    visitor: &mut V,
    current_id: usize,
    visited_nodes: &mut HashSet<usize>,
    visited_edges: &mut HashSet<(usize, usize)>,
) -> Result<(), SmilesError> {
    if visited_nodes.contains(&current_id) {
        return Ok(());
    }

    visited_nodes.insert(current_id);
    visitor.enter_node(smiles, current_id)?;

    let bonds = smiles.edges_for_node(current_id);
    let mut unvisited_bonds = Vec::new();

    for bond in bonds {
        if let Some(other_id) = bond.other(current_id) {
            let edge_key = Smiles::edge_key(current_id, other_id);

            if visited_edges.contains(&edge_key) {
                continue;
            }
            if visited_nodes.contains(&other_id) {
                visited_edges.insert(edge_key);
                visitor.cycle_edge(smiles, current_id, other_id, bond.bond())?;
            } else {
                unvisited_bonds.push((bond, other_id));
            }
        }
    }

    if let Some((main_bond, main_other_id)) = unvisited_bonds.first().copied() {
        visited_edges.insert(Smiles::edge_key(current_id, main_other_id));
        visitor.tree_edge(smiles, *main_bond)?;
        dfs(smiles, visitor, main_other_id, visited_nodes, visited_edges)?;

        for (branch_bond, branch_other_id) in unvisited_bonds.into_iter().skip(1) {
            if visited_nodes.contains(&branch_other_id) {
                continue;
            }
            visited_edges.insert(Smiles::edge_key(current_id, branch_other_id));
            visitor.open_branch(smiles, current_id, branch_other_id)?;
            visitor.tree_edge(smiles, *branch_bond)?;
            dfs(smiles, visitor, branch_other_id, visited_nodes, visited_edges)?;
            visitor.close_branch(smiles, current_id, branch_other_id)?;
        }
    }

    visitor.exit_node(smiles, current_id)?;
    Ok(())
}
