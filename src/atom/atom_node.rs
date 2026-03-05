//! Module for the structure of an atom as a node for use in a `Smiles` graph

use std::{cmp::Ordering, ops::Range};

use crate::{atom::Atom, bond::ring_num::RingNum};
/// Contains information about atom parsed from the SMILES string
#[derive(Clone, PartialEq, Eq)]
pub struct AtomNode {
    /// Unique identifier for each node
    id: usize,
    /// Atom
    atom: Atom,
    /// Span for the atom from the original string
    span: Range<usize>,
    /// Possible Ring Number for this node
    ring_num: Option<RingNum>,
}

impl AtomNode {
    /// Creates a new node
    #[must_use]
    pub fn new(atom: Atom, id: usize, span: Range<usize>, ring_num: Option<RingNum>) -> Self {
        Self { id, atom, span, ring_num }
    }
    /// returns the id
    #[must_use]
    pub fn id(&self) -> usize {
        self.id
    }
    /// returns a borrowed [`Atom`]
    #[must_use]
    pub fn atom(&self) -> &Atom {
        &self.atom
    }
    /// returns the borrowed `Option<RingNum>`
    #[must_use]
    pub fn ring_num(&self) -> &Option<RingNum> {
        &self.ring_num
    }
    /// returns the [`RingNum`] value
    #[must_use]
    pub fn ring_num_val(&self) -> Option<u8> {
        self.ring_num.map(|num| num.get())
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
    use std::{cmp::Ordering, ops::Range};

    use crate::atom::{
        Atom, atom_node::AtomNode, atom_symbol::AtomSymbol, bracketed::BracketAtom,
        unbracketed::UnbracketedAtom,
    };

    #[test]
    fn test_atom_node_all_fields_and_implementations() {
        let ac_symbol = AtomSymbol::new(Some(elements_rs::Element::Ac));
        let atom: Atom = BracketAtom::builder().with_symbol(ac_symbol).build().into();
        let span: Range<usize> = Range { start: 0, end: 1 };
        let atom_node = AtomNode::new(atom.clone(), 0, span, None);
        assert_eq!(&atom_node.id, &0);
        assert_eq!(&atom_node.atom, &atom);
        let next_atom: Atom = UnbracketedAtom::new(AtomSymbol::WildCard, false).into();
        let next_range: Range<usize> = Range { start: 1, end: 2 };
        let next_atom_node = AtomNode::new(next_atom.clone(), 1, next_range, None);
        assert_eq!(atom_node.cmp(&next_atom_node), Ordering::Less);
    }
}
