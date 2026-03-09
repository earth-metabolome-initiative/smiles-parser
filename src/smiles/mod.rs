//! Represents a SMILES structure.

use std::{collections::BTreeSet, fmt};

use crate::{
    atom::atom_node::AtomNode,
    bond::{Bond, bond_edge::BondEdge},
    errors::SmilesError,
};

mod from_str;

/// Represents a SMILES structure.
pub struct Smiles {
    atom_nodes: Vec<AtomNode>,
    bond_edges: Vec<BondEdge>,
}

impl Smiles {
    /// creates a new instance
    #[must_use]
    pub fn new() -> Self {
        Self { atom_nodes: Vec::new(), bond_edges: Vec::new() }
    }
    /// Pushes an AtomNode
    pub fn push_node(&mut self, node: AtomNode) {
        self.atom_nodes.push(node);
    }
    /// adds an edge from two nodes and the [`Bond`]
    ///
    /// # Errors
    /// - Returns a [`SmilesError::NodeIdInvalid`] if a node cannot be found in
    ///   the edge list
    pub fn push_edge(
        &mut self,
        node_a: usize,
        node_b: usize,
        bond: Bond,
    ) -> Result<(), SmilesError> {
        self.atom_nodes.sort();
        // use the NodeIdInvalid for err
        if self
            .atom_nodes
            .binary_search_by_key(&node_a, super::atom::atom_node::AtomNode::id)
            .is_err()
        {
            return Err(SmilesError::NodeIdInvalid(node_a));
        }
        if self
            .atom_nodes
            .binary_search_by_key(&node_b, super::atom::atom_node::AtomNode::id)
            .is_err()
        {
            return Err(SmilesError::NodeIdInvalid(node_b));
        }
        let bond_edge = BondEdge::new(node_a, node_b, bond);
        self.bond_edges.push(bond_edge);

        Ok(())
    }
    /// Returns slice of the nodes
    #[must_use]
    pub fn nodes(&self) -> &[AtomNode] {
        &self.atom_nodes
    }
    /// Returns mutable slice of nodes
    #[must_use]
    pub fn nodes_mut(&mut self) -> &mut [AtomNode] {
        &mut self.atom_nodes
    }
    /// Returns a node if exists from an `id`
    #[must_use]
    pub fn node_by_id(&self, id: usize) -> Option<&AtomNode> {
        self.atom_nodes.iter().find(|node| node.id() == id)
    }
    /// Returns slice of the edges
    #[must_use]
    pub fn edges(&self) -> &[BondEdge] {
        &self.bond_edges
    }
    /// Returns mutable slice of edges
    pub fn edges_mut(&mut self) -> &mut [BondEdge] {
        &mut self.bond_edges
    }
    /// Find all neighboring edges bonded to the current node id
    #[must_use]
    pub fn neighbors(&self, id: usize) -> Vec<(usize, Bond)> {
        self.bond_edges
            .iter()
            .filter_map(|edge| {
                if edge.node_a() == id {
                    Some((edge.node_b(), edge.bond().to_owned()))
                } else if edge.node_b() == id {
                    Some((edge.node_a(), edge.bond().to_owned()))
                } else {
                    None
                }
            })
            .collect()
    }
}

impl Default for Smiles {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Smiles {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fn write_subtree(
            smiles: &Smiles,
            f: &mut fmt::Formatter<'_>,
            current: usize,
            parent: Option<usize>,
            visited: &mut BTreeSet<usize>,
        ) -> fmt::Result {
            visited.insert(current);

            let node = smiles.node_by_id(current).ok_or(fmt::Error)?;
            write!(f, "{node}")?;

            let mut children = smiles
                .neighbors(current)
                .into_iter()
                .filter(|(next, _)| Some(*next) != parent && !visited.contains(next))
                .collect::<Vec<_>>();

            children.sort_by_key(|(id, _)| {
                smiles.node_by_id(*id).map_or(usize::MAX, |node| node.span().start)
            });

            if let Some((first_child, first_bond)) = children.first().copied() {
                if !matches!(first_bond, Bond::Single) {
                    write!(f, "{first_bond}")?;
                }
                write_subtree(smiles, f, first_child, Some(current), visited)?;

                for (child, bond) in children.into_iter().skip(1) {
                    f.write_str("(")?;
                    if !matches!(bond, Bond::Single) {
                        write!(f, "{bond}")?;
                    }
                    write_subtree(smiles, f, child, Some(current), visited)?;
                    f.write_str(")")?;
                }
            }

            Ok(())
        }

        let mut roots = self.atom_nodes.iter().collect::<Vec<_>>();
        roots.sort_by_key(|node| node.span().start);

        let mut visited = BTreeSet::new();
        let mut first_component = true;

        for node in roots {
            if visited.contains(&node.id()) {
                continue;
            }

            if !first_component {
                f.write_str(".")?;
            }
            first_component = false;

            write_subtree(self, f, node.id(), None, &mut visited)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use elements_rs::Element;

    use crate::{
        atom::{Atom, atom_node::AtomNode, atom_symbol::AtomSymbol, unbracketed::UnbracketedAtom},
        bond::Bond,
        errors::SmilesError,
        smiles::Smiles,
    };

    fn atom_node(id: usize, element: Element) -> AtomNode {
        let atom: Atom = UnbracketedAtom::new(AtomSymbol::Element(element), false).into();
        AtomNode::new(atom, id, id..(id + 1), None)
    }

    #[test]
    fn test_smiles_display_single_atom() {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C));

        assert_eq!(smiles.to_string(), "C");
    }

    #[test]
    fn test_smiles_display_linear_single_bonds_are_implicit() -> Result<(), SmilesError> {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C));
        smiles.push_node(atom_node(1, Element::C));
        smiles.push_node(atom_node(2, Element::O));

        smiles.push_edge(0, 1, Bond::Single)?;
        smiles.push_edge(1, 2, Bond::Single)?;

        assert_eq!(smiles.to_string(), "CCO");
        Ok(())
    }

    #[test]
    fn test_smiles_display_explicit_multiple_bond() -> Result<(), SmilesError> {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C));
        smiles.push_node(atom_node(1, Element::O));

        smiles.push_edge(0, 1, Bond::Double)?;

        assert_eq!(smiles.to_string(), "C=O");
        Ok(())
    }

    #[test]
    fn test_smiles_display_branching() -> Result<(), SmilesError> {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C));
        smiles.push_node(atom_node(1, Element::O));
        smiles.push_node(atom_node(2, Element::N));

        smiles.push_edge(0, 1, Bond::Single)?;
        smiles.push_edge(0, 2, Bond::Single)?;

        assert_eq!(smiles.to_string(), "CO(N)");
        Ok(())
    }

    #[test]
    fn test_smiles_display_branch_with_explicit_bond() -> Result<(), SmilesError> {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C));
        smiles.push_node(atom_node(1, Element::O));
        smiles.push_node(atom_node(2, Element::N));

        smiles.push_edge(0, 1, Bond::Double)?;
        smiles.push_edge(0, 2, Bond::Single)?;

        assert_eq!(smiles.to_string(), "C=O(N)");
        Ok(())
    }

    #[test]
    fn test_smiles_display_disconnected_components() {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C));
        smiles.push_node(atom_node(1, Element::O));

        assert_eq!(smiles.to_string(), "C.O");
    }

    #[test]
    fn test_smiles_display_branch_off_middle_of_chain() -> Result<(), SmilesError> {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C));
        smiles.push_node(atom_node(1, Element::C));
        smiles.push_node(atom_node(2, Element::O));
        smiles.push_node(atom_node(3, Element::N));

        smiles.push_edge(0, 1, Bond::Single)?;
        smiles.push_edge(1, 2, Bond::Single)?;
        smiles.push_edge(1, 3, Bond::Single)?;

        assert_eq!(smiles.to_string(), "CCO(N)");
        Ok(())
    }
}
