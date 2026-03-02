//! Module for the structure of an atom as a node for use in a [`Smiles`] graph

use std::cmp::Ordering;

use crate::atom::Atom;
/// Contains information about atom parsed from the SMILES string
#[derive(PartialEq, Eq)]
pub struct AtomNode {
    /// Unique identifier for each node
    id: usize,
    /// Atom
    atom: Atom,
}

impl AtomNode {
    /// Creates a new node
    #[must_use]
    pub fn new(atom: Atom, id: usize) -> Self {
        Self { id, atom }
    }
    /// returns the id
    #[must_use]
    pub fn id(&self) -> usize {
        self.id
    }
    /// returns the [`Atom`]
    #[must_use]
    pub fn atom(&self) -> &Atom {
        &self.atom
    }
}

impl PartialOrd for AtomNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AtomNode {
    fn cmp(&self, other: &Self) -> Ordering {
        self.id.cmp(&other.id)
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use crate::atom::{Atom, atom_node::AtomNode, atom_symbol::AtomSymbol, bracketed::BracketAtom, unbracketed::UnbracketedAtom};

    #[test]
    fn test_atom_node_all_fields_and_implementations() {
        let ac_symbol = AtomSymbol::new(Some(elements_rs::Element::Ac));
        let atom: Atom = BracketAtom::builder().with_symbol(ac_symbol).build().into();
        let atom_node = AtomNode::new(atom.clone(), 0);
        assert_eq!(&atom_node.id, &0);
        assert_eq!(&atom_node.atom, &atom);
        let next_atom: Atom = UnbracketedAtom::new(AtomSymbol::WildCard, false).into();
        let next_atom_node = AtomNode::new(next_atom.clone(), 1);
        assert_eq!(atom_node.cmp(&next_atom_node), Ordering::Less);

    }
}