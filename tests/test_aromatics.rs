//! Tests on Elements that should or should not be parsed as aromatic

use smiles_parser::{parser::token_iter::TokenIter, token::Token};
const SMILES_STR: &[&str] = &[
    "c1ccccc1",          // benzene
    "n1ccccc1",          // pyridine
    "Cc1ccccc1",         // toluene
    "c1cccc2ccccc12",    // naphthalene
    "n1cc[nH]c1",        // imidazole (bracketed aromatic N)
    "C1CCCCC1",          // cyclohexane
    "C1CCC=CC1",         // cyclohexene
    "C1CCC1",            // cyclobutane
    "CCO",               // ethanol
    "C1CCCC2CCCCC12",    // decalin
    "[nH]1cccc1",        // pyrrole-like aromatic N
    "c1[n+]([O-])cccc1", // nitro-substituted aromatic ring
    "[CH3][CH2][OH]",    // ethanol (explicit atoms)
    "[NH4+]",            // ammonium (no aromaticity)
];

#[test]
fn test_aromatic_benzene_from_tokenization() {
    //  c1ccccc1
    let benzene_tokens = vec![
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::Label(1),
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::Label(1),
    ];
    let benzene_line = SMILES_STR[0];
    let tokens = TokenIter::from(benzene_line)
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_else(|_| panic!("Failed to parse {benzene_line}"));
    assert_eq!(benzene_tokens, tokens);
}

#[test]
fn test_aromatic_imidazole_from_tokenization() {
    //  n1cc[nH]c1
    let imidazole_tokens = vec![
        Token::Atom { element: elements_rs::Element::N, aromatic: true },
        Token::Label(1),
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::LeftSquareBracket,
        Token::Atom { element: elements_rs::Element::N, aromatic: true },
        Token::Atom { element: elements_rs::Element::H, aromatic: false },
        Token::RightSquareBracket,
        Token::Atom { element: elements_rs::Element::C, aromatic: true },
        Token::Label(1),
    ];
    let imidazole_line = SMILES_STR[4];
    let tokens = TokenIter::from(imidazole_line)
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_else(|_| panic!("Failed to parse {imidazole_line}"));
    assert_eq!(imidazole_tokens, tokens);
}
