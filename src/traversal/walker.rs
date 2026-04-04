//! Module for walking the [`Smiles`] graph and tracking progress

use alloc::{
    collections::{BTreeMap, BTreeSet},
    vec::Vec,
};

use crate::{
    bond::bond_edge::BondEdge, errors::SmilesError, smiles::Smiles,
    traversal::visitor_trait::Visitor,
};

/// Traverses the [`Smiles`] graph via a depth first search
///
/// # Errors
/// - Will return [`SmilesError::NodeIdInvalid`] if a node is unable to be found
pub fn walk<V: Visitor>(smiles: &Smiles, visitor: &mut V) -> Result<(), SmilesError> {
    // Each key is a node ID and value is the list of incident edge indices.
    let mut adjacency: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    build_adjacency(&mut adjacency, smiles.edges());
    let mut visited_nodes: BTreeSet<usize> = BTreeSet::new();
    let mut visited_edges: BTreeSet<usize> = BTreeSet::new();
    let mut component_index = 0;

    for node in smiles.nodes() {
        if !visited_nodes.contains(&node.id()) {
            visitor.start_component(smiles, node.id(), component_index)?;
            dfs(smiles, visitor, node.id(), &mut visited_nodes, &mut visited_edges, &adjacency)?;
            visitor.finish_component(smiles, node.id(), component_index)?;
            component_index += 1;
        }
    }

    Ok(())
}

fn build_adjacency(adjacency: &mut BTreeMap<usize, Vec<usize>>, edges: &[BondEdge]) {
    for (index, edge) in edges.iter().enumerate() {
        let node_a = edge.node_a();
        let node_b = edge.node_b();
        adjacency.entry(node_a).or_default().push(index);
        adjacency.entry(node_b).or_default().push(index);
    }
}

fn dfs<V: Visitor>(
    smiles: &Smiles,
    visitor: &mut V,
    current_id: usize,
    visited_nodes: &mut BTreeSet<usize>,
    visited_edges: &mut BTreeSet<usize>,
    adjacency: &BTreeMap<usize, Vec<usize>>,
) -> Result<(), SmilesError> {
    if !visited_nodes.insert(current_id) {
        return Ok(());
    }

    visitor.enter_node(smiles, current_id)?;

    let mut unvisited_bonds = Vec::new();
    let mut queued_neighbors = BTreeSet::new();

    if let Some(bond_indices) = adjacency.get(&current_id) {
        for bond_index in bond_indices {
            if visited_edges.contains(bond_index) {
                continue;
            }
            let bond = &smiles.edges()[*bond_index];
            if let Some(other_id) = bond.other(current_id) {
                if visited_nodes.contains(&other_id) {
                    visited_edges.insert(*bond_index);
                    visitor.cycle_edge(smiles, current_id, other_id, bond.bond())?;
                } else if queued_neighbors.insert(other_id) {
                    unvisited_bonds.push((*bond_index, bond, other_id));
                }
            }
        }
    }

    // Prefer a leaf-like continuation edge when multiple choices exist. This
    // keeps render output stable for graphs where one root-adjacent edge is
    // also reachable through a later cycle closure.
    unvisited_bonds.sort_by_key(|(_, _, other_id)| {
        let degree = adjacency.get(other_id).map_or(0, Vec::len);
        (degree == 1, *other_id)
    });

    if let Some((main_edge_index, main_bond, main_other_id)) = unvisited_bonds.pop() {
        for (branch_edge_index, branch_bond, branch_other_id) in unvisited_bonds {
            if visited_nodes.contains(&branch_other_id) {
                continue;
            }
            visited_edges.insert(branch_edge_index);
            visitor.open_branch(smiles, current_id, branch_other_id)?;
            visitor.tree_edge(smiles, *branch_bond)?;
            dfs(smiles, visitor, branch_other_id, visited_nodes, visited_edges, adjacency)?;
            visitor.close_branch(smiles, current_id, branch_other_id)?;
        }
        if !visited_nodes.contains(&main_other_id) {
            visited_edges.insert(main_edge_index);
            visitor.tree_edge(smiles, *main_bond)?;
            dfs(smiles, visitor, main_other_id, visited_nodes, visited_edges, adjacency)?;
        }
    }
    visitor.exit_node(smiles, current_id)?;
    Ok(())
}
