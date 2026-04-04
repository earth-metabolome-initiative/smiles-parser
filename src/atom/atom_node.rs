//! Module for the structure of an atom as a node for use in a `Smiles` graph

use core::{fmt, ops::Range};

use crate::atom::Atom;
/// Contains information about atom parsed from the SMILES string
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomNode {
    /// Unique identifier for each node
    id: usize,
    /// Atom
    atom: Atom,
    /// Span for the atom from the original string
    span: Range<usize>,
}

impl AtomNode {
    /// Creates a new node
    #[must_use]
    pub fn new(atom: Atom, id: usize, span: Range<usize>) -> Self {
        Self { id, atom, span }
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
    /// returns the `Span` of the node
    #[must_use]
    pub fn span(&self) -> &Range<usize> {
        &self.span
    }
}

impl fmt::Display for AtomNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.atom)
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;
    use core::ops::Range;

    use elements_rs::Element;

    use crate::atom::{
        Atom,
        atom_node::AtomNode,
        atom_symbol::AtomSymbol,
        bracketed::{BracketAtom, charge::Charge, hydrogen_count::HydrogenCount},
        unbracketed::UnbracketedAtom,
    };

    #[test]
    fn test_atom_node_new_and_accessors() {
        let atom: Atom =
            BracketAtom::builder().with_symbol(AtomSymbol::Element(Element::C)).build().into();
        let span = Range { start: 2, end: 5 };

        let node = AtomNode::new(atom.clone(), 42, span.clone());

        assert_eq!(node.id(), 42);
        assert_eq!(node.atom(), &atom);
        assert_eq!(node.span(), &span);
    }

    #[test]
    fn test_atom_node_ring_num_none() {
        let atom: Atom = UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into();
        let node = AtomNode::new(atom, 0, 0..1);

        assert_eq!(node.span(), &(0..1));
    }

    #[test]
    fn test_atom_node_display_delegates_to_atom_unbracketed() {
        let atom: Atom = UnbracketedAtom::new(AtomSymbol::Element(Element::C), true).into();
        let node = AtomNode::new(atom, 0, 0..1);

        assert_eq!(node.to_string(), "c");
    }

    #[test]
    fn test_atom_node_display_delegates_to_atom_bracketed() {
        let atom: Atom = BracketAtom::builder()
            .with_symbol(AtomSymbol::Element(Element::C))
            .with_isotope(13)
            .with_hydrogens(HydrogenCount::Explicit(2))
            .with_charge(Charge::try_new(-1).unwrap())
            .build()
            .into();

        let node = AtomNode::new(atom, 1, 3..9);

        assert_eq!(node.to_string(), "[13CH2-]");
    }

    #[test]
    fn test_atom_node_equality_uses_all_fields() {
        let atom: Atom = UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into();

        let a = AtomNode::new(atom.clone(), 1, 0..1);
        let b = AtomNode::new(atom.clone(), 1, 0..1);
        let c = AtomNode::new(atom, 1, 1..2);

        assert_eq!(a, b);
        assert_ne!(a, c);
    }
}
