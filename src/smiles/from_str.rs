use core::str::FromStr;

use super::Smiles;
use crate::{errors::SmilesErrorWithSpan, parser::smiles_parser::parse_smiles};

impl FromStr for Smiles {
    type Err = SmilesErrorWithSpan;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_smiles(s)
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;
    use std::str::FromStr;

    use elements_rs::Element;

    use crate::{
        atom::{Atom, atom_symbol::AtomSymbol},
        bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
        smiles::Smiles,
    };

    #[test]
    fn parse_benzene_with_ring_nums() {
        let smiles = Smiles::from_str("C1=CC=CC=C1")
            .unwrap_or_else(|e| panic!("Failed to tokenize:\n{}", e.render("C1=CC=CC=C1")));

        let expected_nodes = [
            Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
            Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
            Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
            Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
            Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
            Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
        ];
        let expected_ring_edge =
            BondEdge::new(0, 5, Bond::Single, Some(RingNum::try_new(1).unwrap()));

        assert_eq!(smiles.nodes()[0], expected_nodes[0]);
        assert_eq!(smiles.nodes()[5], expected_nodes[5]);

        assert_eq!(smiles.edge_for_node_pair((0, 5)), Some(expected_ring_edge));
        assert_eq!(smiles.edge_for_node_pair((0, 5)).and_then(|edge| edge.ring_num_val()), Some(1));
        assert_eq!(smiles.node_by_id(0), Some(&expected_nodes[0]));
        assert_eq!(smiles.node_by_id(5), Some(&expected_nodes[5]));
    }

    #[test]
    fn from_str_propagates_token_iter_error() {
        let err = Smiles::from_str("Ac").expect_err("expected tokenization to fail");

        assert_eq!(
            err.smiles_error(),
            crate::errors::SmilesError::InvalidUnbracketedAtom(AtomSymbol::Element(Element::Ac))
        );
        assert_eq!(err.start(), 0);
        assert_eq!(err.end(), 2);
        assert_eq!(err.to_string(), "Invalid unbracketed atom: Ac at 0..2");
    }

    #[test]
    fn from_str_propagates_parser_error() {
        let err = Smiles::from_str("C(").expect_err("expected parsing to fail");

        assert_eq!(err.smiles_error(), crate::errors::SmilesError::UnclosedBranch);
        assert_eq!(err.start(), 1);
        assert_eq!(err.end(), 2);
        assert_eq!(err.to_string(), "Branch not closed at 1..2");
    }
}
