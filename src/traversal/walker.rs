//! Module for walking the [`Smiles`] graph and tracking progress.

use alloc::vec::Vec;

use geometric_traits::traits::{SizedRowsSparseMatrix2D, SparseMatrix2D, SparseValuedMatrix2DRef};

use crate::{errors::SmilesError, smiles::Smiles, traversal::visitor_trait::Visitor};

/// Traverses the [`Smiles`] graph via a depth first search.
///
/// # Errors
/// - Will return [`SmilesError::NodeIdInvalid`] if a node is unable to be
///   found.
pub(crate) fn walk<V: Visitor>(smiles: &Smiles, visitor: &mut V) -> Result<(), SmilesError> {
    let node_count = smiles.nodes().len();
    let mut visited_nodes = vec![false; node_count];
    let mut visited_edges = vec![false; smiles.number_of_bonds()];
    let mut component_index = 0;

    for node_id in 0..node_count {
        if !visited_nodes[node_id] {
            visitor.start_component(smiles, node_id, component_index)?;
            dfs(smiles, visitor, node_id, &mut visited_nodes, &mut visited_edges)?;
            visitor.finish_component(smiles, node_id, component_index)?;
            component_index += 1;
        }
    }

    Ok(())
}

fn dfs<V: Visitor>(
    smiles: &Smiles,
    visitor: &mut V,
    current_id: usize,
    visited_nodes: &mut [bool],
    visited_edges: &mut [bool],
) -> Result<(), SmilesError> {
    if visited_nodes[current_id] {
        return Ok(());
    }
    visited_nodes[current_id] = true;

    visitor.enter_node(smiles, current_id)?;

    let bond_matrix = smiles.bond_matrix();
    let mut unvisited_bonds =
        Vec::with_capacity(bond_matrix.number_of_defined_values_in_row(current_id));

    for (other_id, entry) in
        bond_matrix.sparse_row(current_id).zip(bond_matrix.sparse_row_values_ref(current_id))
    {
        let edge_order = entry.order();
        if visited_edges[edge_order] {
            continue;
        }

        if visited_nodes[other_id] {
            visited_edges[edge_order] = true;
            visitor.cycle_edge(smiles, current_id, other_id, entry.bond())?;
        } else {
            let is_leaf = bond_matrix.number_of_defined_values_in_row(other_id) == 1;
            unvisited_bonds.push((other_id, is_leaf, *entry));
        }
    }

    // Prefer a leaf-like continuation edge when multiple choices exist. This
    // keeps render output stable for graphs where one root-adjacent edge is
    // also reachable through a later cycle closure.
    unvisited_bonds.sort_by_key(|(other_id, is_leaf, entry)| (*is_leaf, entry.order(), *other_id));

    if let Some((main_other_id, _, main_entry)) = unvisited_bonds.pop() {
        for (branch_other_id, _, branch_entry) in unvisited_bonds {
            if visited_nodes[branch_other_id] {
                continue;
            }
            visited_edges[branch_entry.order()] = true;
            visitor.open_branch(smiles, current_id, branch_other_id)?;
            visitor.tree_edge(smiles, branch_entry.to_bond_edge(current_id, branch_other_id))?;
            dfs(smiles, visitor, branch_other_id, visited_nodes, visited_edges)?;
            visitor.close_branch(smiles, current_id, branch_other_id)?;
        }
        if !visited_nodes[main_other_id] {
            visited_edges[main_entry.order()] = true;
            visitor.tree_edge(smiles, main_entry.to_bond_edge(current_id, main_other_id))?;
            dfs(smiles, visitor, main_other_id, visited_nodes, visited_edges)?;
        }
    }

    visitor.exit_node(smiles, current_id)?;
    Ok(())
}
