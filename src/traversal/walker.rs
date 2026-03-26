//! Module for walking the [`Smiles`] graph and tracking progress

use std::collections::HashSet;

use crate::{errors::SmilesError, smiles::Smiles, traversal::visitor_trait::Visitor};

/// Traverses the [`Smiles`] graph via a depth first search
///
/// # Errors
/// - Will return [`SmilesError::NodeIdInvalid`] if a node is unable to be found
pub fn walk<V: Visitor>(smiles: &Smiles, visitor: &mut V) -> Result<(), SmilesError> {
    let mut visited_nodes: HashSet<usize> = HashSet::new();
    let mut visited_edges: HashSet<usize> = HashSet::new();
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
    visited_edges: &mut HashSet<usize>,
) -> Result<(), SmilesError> {
    if visited_nodes.contains(&current_id) {
        return Ok(());
    }

    visited_nodes.insert(current_id);
    visitor.enter_node(smiles, current_id)?;

    let mut unvisited_bonds = Vec::new();
    let mut queued_neighbors = HashSet::new();

    for (edge_index, bond) in smiles.edges().iter().enumerate() {
        if !bond.contains(current_id) || visited_edges.contains(&edge_index) {
            continue;
        }
        if let Some(other_id) = bond.other(current_id) {
            if visited_nodes.contains(&other_id) {
                visited_edges.insert(edge_index);
                visitor.cycle_edge(smiles, current_id, other_id, bond.bond())?;
            } else if queued_neighbors.insert(other_id) {
                unvisited_bonds.push((edge_index, bond, other_id));
            }
        }
    }

    if let Some((main_edge_index, main_bond, main_other_id)) = unvisited_bonds.pop() {
        for (branch_edge_index, branch_bond, branch_other_id) in unvisited_bonds {
            if visited_nodes.contains(&branch_other_id) {
                continue;
            }
            visited_edges.insert(branch_edge_index);
            visitor.open_branch(smiles, current_id, branch_other_id)?;
            visitor.tree_edge(smiles, *branch_bond)?;
            dfs(smiles, visitor, branch_other_id, visited_nodes, visited_edges)?;
            visitor.close_branch(smiles, current_id, branch_other_id)?;
        }
        if !visited_nodes.contains(&main_other_id) {
            visited_edges.insert(main_edge_index);
            visitor.tree_edge(smiles, *main_bond)?;
            dfs(smiles, visitor, main_other_id, visited_nodes, visited_edges)?;
        }
    }
    visitor.exit_node(smiles, current_id)?;
    Ok(())
}
