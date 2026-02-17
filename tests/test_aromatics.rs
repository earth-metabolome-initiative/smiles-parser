//! Tests on Elements that should or should not be parsed as aromatic

use elements_rs::Element;
use smiles_parser::{
    atom_symbol::AtomSymbol, errors::SmilesError, parser::token_iter::TokenIter, ring_num::RingNum,
    token::Token, unbracketed::UnbracketedAtom,
};
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
fn test_aromatic_benzene_from_tokenization() -> Result<(), SmilesError> {
    // c1ccccc1
    let aromatic_c = UnbracketedAtom::new(AtomSymbol::Element(Element::C), true);

    let expected = vec![
        Token::UnbracketedAtom(aromatic_c),
        Token::RingClosure(RingNum::try_new(1)?),
        Token::UnbracketedAtom(aromatic_c),
        Token::UnbracketedAtom(aromatic_c),
        Token::UnbracketedAtom(aromatic_c),
        Token::UnbracketedAtom(aromatic_c),
        Token::UnbracketedAtom(aromatic_c),
        Token::RingClosure(RingNum::try_new(1)?),
    ];

    let line = SMILES_STR[0];
    let got = TokenIter::from(line).collect::<Result<Vec<_>, _>>()?;
    assert_eq!(expected, got);

    Ok(())
}

#[test]
fn test_aromatic_imidazole_from_tokenization() -> Result<(), SmilesError> {
    // n1cc[nH]c1
    let aromatic_n = UnbracketedAtom::new(AtomSymbol::Element(Element::N), true);
    let aromatic_c = UnbracketedAtom::new(AtomSymbol::Element(Element::C), true);

    let expected = vec![
        Token::UnbracketedAtom(aromatic_n),
        Token::RingClosure(RingNum::try_new(1)?),
        Token::UnbracketedAtom(aromatic_c),
        Token::UnbracketedAtom(aromatic_c),
        Token::BracketedAtom(
            smiles_parser::bracketed::bracket_atom::BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::N))
                .with_aromatic(true)
                .with_hydrogens(smiles_parser::bracketed::hydrogen_count::HydrogenCount::new(Some(
                    1,
                )))
                .build(),
        ),
        Token::UnbracketedAtom(aromatic_c),
        Token::RingClosure(RingNum::try_new(1)?),
    ];

    let line = SMILES_STR[4];
    let got = TokenIter::from(line).collect::<Result<Vec<_>, _>>()?;
    assert_eq!(expected, got);
    Ok(())
}
