use std::str::FromStr;

use super::Smiles;
use crate::{
    errors::SmilesErrorWithSpan,
    parser::{smiles_parser::SmilesParser, token_iter::TokenIter},
};

impl FromStr for Smiles {
    type Err = SmilesErrorWithSpan;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let token_iter = TokenIter::from(s);
        let tokens = token_iter.collect::<Result<Vec<_>, _>>()?;
        SmilesParser::new(&tokens).parse()
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use elements_rs::Element;

    use crate::{
        atom::{atom_node::AtomNode, atom_symbol::AtomSymbol, unbracketed::UnbracketedAtom},
        bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
        smiles::Smiles,
    };

    #[test]
    fn parse_benzene_with_ring_nums() {
        let smiles = Smiles::from_str("C1=CC=CC=C1")
            .unwrap_or_else(|e| panic!("Failed to tokenize:\n{}", e.render("C1=CC=CC=C1")));

        let expected = Smiles {
            atom_nodes: vec![
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    0,
                    0..1,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    1,
                    3..4,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    2,
                    4..5,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    3,
                    6..7,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    4,
                    7..8,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    5,
                    9..10,
                ),
            ],
            bond_edges: vec![
                BondEdge::new(0, 1, Bond::Double, None),
                BondEdge::new(1, 2, Bond::Single, None),
                BondEdge::new(2, 3, Bond::Double, None),
                BondEdge::new(3, 4, Bond::Single, None),
                BondEdge::new(4, 5, Bond::Double, None),
                BondEdge::new(5, 0, Bond::Single, Some(RingNum::try_new(1).unwrap())),
            ],
        };

        assert_eq!(smiles.nodes()[0], expected.atom_nodes[0]);
        assert_eq!(smiles.nodes()[5], expected.atom_nodes[5]);

        assert_eq!(smiles.edges()[5], expected.bond_edges[5]);
        assert_eq!(smiles.edges()[5].ring_num_val(), Some(1));
        assert_eq!(smiles.node_by_id(0), expected.node_by_id(0));
        assert_eq!(smiles.node_by_id(5), expected.node_by_id(5));
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
