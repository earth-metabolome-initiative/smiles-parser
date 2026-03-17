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
                    Some(RingNum::try_new(1).unwrap()),
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    1,
                    3..4,
                    None,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    2,
                    4..5,
                    None,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    3,
                    6..7,
                    None,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    4,
                    7..8,
                    None,
                ),
                AtomNode::new(
                    UnbracketedAtom::new(AtomSymbol::Element(Element::C), false).into(),
                    5,
                    9..10,
                    None,
                ),
            ],
            bond_edges: vec![
                BondEdge::new(0, 1, Bond::Double),
                BondEdge::new(1, 2, Bond::Single),
                BondEdge::new(2, 3, Bond::Double),
                BondEdge::new(3, 4, Bond::Single),
                BondEdge::new(4, 5, Bond::Double),
                BondEdge::new(5, 0, Bond::Single),
            ],
        };
        assert_eq!(smiles.atom_nodes[0], expected.atom_nodes[0]);
        assert_eq!(smiles.atom_nodes[0].ring_num_val(), expected.atom_nodes[0].ring_num_val());
        assert_eq!(smiles.atom_nodes[5], expected.atom_nodes[5]);
        assert_eq!(smiles.atom_nodes[5].ring_num_val(), Some(1));
    }
}
