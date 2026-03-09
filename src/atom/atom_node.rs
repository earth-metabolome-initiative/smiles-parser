//! Module for the structure of an atom as a node for use in a `Smiles` graph

use std::{cmp::Ordering, fmt, ops::Range};

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
    /// returns the [`Span`] of the node
    #[must_use]
    pub fn span(&self) -> &Range<usize> {
        &self.span
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

impl fmt::Display for AtomNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.atom)
    }
}

#[cfg(test)]
mod tests {
    use std::{cmp::Ordering, ops::Range};

    use elements_rs::Element;

    use crate::{
        atom::{
            Atom,
            atom_node::AtomNode,
            atom_symbol::AtomSymbol,
            bracketed::{BracketAtom, charge::Charge, hydrogen_count::HydrogenCount},
            unbracketed::UnbracketedAtom,
        },
        bond::ring_num::RingNum,
    };

    #[test]
    fn test_atom_node_new_and_accessors() {
        let atom: Atom =
            BracketAtom::builder().with_symbol(AtomSymbol::Element(Element::C)).build().into();
        let span = Range { start: 2, end: 5 };
        let ring_num = Some(RingNum::try_new(7).unwrap());

        let node = AtomNode::new(atom.clone(), 42, span.clone(), ring_num);

        assert_eq!(node.id(), 42);
        assert_eq!(node.atom(), &atom);
        assert_eq!(node.span(), &span);
        assert_eq!(node.ring_num(), &Some(RingNum::try_new(7).unwrap()));
        assert_eq!(node.ring_num_val(), Some(7));
    }

    #[test]
    fn test_atom_node_ring_num_none() {
        let atom: Atom = UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into();
        let node = AtomNode::new(atom, 0, 0..1, None);

        assert_eq!(node.span(), &(0..1));
        assert_eq!(node.ring_num(), &None);
        assert_eq!(node.ring_num_val(), None);
    }

    #[test]
    fn test_atom_node_ord_and_partial_ord_compare_by_id_only() {
        let atom_a: Atom = UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into();
        let atom_b: Atom = BracketAtom::builder()
            .with_symbol(AtomSymbol::Element(Element::O))
            .with_hydrogens(HydrogenCount::Explicit(1))
            .with_charge(Charge::try_new(-1).unwrap())
            .build()
            .into();

        let node_a = AtomNode::new(atom_a, 1, 0..1, None);
        let node_b = AtomNode::new(atom_b, 2, 10..15, Some(RingNum::try_new(12).unwrap()));

        assert_eq!(node_a.cmp(&node_b), Ordering::Less);
        assert_eq!(node_b.cmp(&node_a), Ordering::Greater);
        assert_eq!(node_a.partial_cmp(&node_b), Some(Ordering::Less));
        assert_eq!(node_b.partial_cmp(&node_a), Some(Ordering::Greater));
    }

    #[test]
    fn test_atom_node_ord_equal_when_ids_equal() {
        let atom_a: Atom = UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into();
        let atom_b: Atom = UnbracketedAtom::new(AtomSymbol::WildCard, false).into();

        let node_a = AtomNode::new(atom_a, 5, 0..1, None);
        let node_b = AtomNode::new(atom_b, 5, 99..100, Some(RingNum::try_new(3).unwrap()));

        assert_eq!(node_a.cmp(&node_b), Ordering::Equal);
        assert_eq!(node_a.partial_cmp(&node_b), Some(Ordering::Equal));
    }

    #[test]
    fn test_atom_node_display_delegates_to_atom_unbracketed() {
        let atom: Atom = UnbracketedAtom::new(AtomSymbol::Element(Element::C), true).into();
        let node = AtomNode::new(atom, 0, 0..1, None);

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

        let node = AtomNode::new(atom, 1, 3..9, Some(RingNum::try_new(1).unwrap()));

        assert_eq!(node.to_string(), "[13CH2-]");
    }
}
