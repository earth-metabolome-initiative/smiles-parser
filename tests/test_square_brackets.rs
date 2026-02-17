//! Test for tokenizing square brackets

use smiles_parser::{errors::SmilesError, parser::token_iter::TokenIter, token::Token};

/// const for testing square brackets
const SMILES_WITH_BRACKETS: &[&str] = &[
    "[OH2]",
    "[OH3+]",
    "[Ti+4]",
    "[Co+3]",
    "[Ga+]$[As-]",
    "[Na+].[Cl-]",
    "N[C@@H](C)C(=O)O",
    "N[CH](C)C(=O)O",
    "[14cH]1ccccc1",
    "[14c@H]1ccccc1",
    "[2H]C(Cl)(Cl)Cl",
    "[C@@H](C)(N)C(=O)O",
    "C[C@H](N)C(=O)O",
    "OC(=O)[C@@H](N)C",
    "[K+].C=C.Cl[Pt-](Cl)Cl.O",
    "[Ti++++]",
    "CCN1C[C@]2(COC)CC[C@H](O)[C@@]34[C@@H]5C[C@H]6[C@H](OC)[C@@H]5[C@](O)(C[C@@H]6OC)[C@@](O)([C@@H](OC)[C@H]23)[C@@H]14",
    "C[C@@H]1C[C@@]2(O[C@H]2C)C(=O)O[C@@H]2CCN(C)C/C=C(/COC(=O)[C@]1(C)O)C2=O",
    "CC1=C[C@H](O)CC(C)(C)[C@H]1/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(\\C)CO",
];
#[test]
fn test_valid_square_brackets() {
    for &s in SMILES_WITH_BRACKETS {
        let _tokens = TokenIter::from(s)
            .collect::<Result<Vec<_>, _>>()
            .unwrap_or_else(|_| panic!("Failed to parse {s}"));
    }
}

#[test]
fn test_unexpected_right_brackets() {
    // New tokenizer does not have a dedicated "UnexpectedRightBracket" error:
    // stray ']' becomes an UnexpectedCharacter.
    let error = TokenIter::from("[Co+3]]").collect::<Result<Vec<_>, SmilesError>>();
    assert!(error.is_err());
    assert_eq!(error, Err(SmilesError::UnexpectedCharacter { character: ']' }));
}

#[test]
fn test_unexpected_left_brackets() {
    let error = TokenIter::from("[[Co+3]").collect::<Result<Vec<_>, SmilesError>>();
    assert!(error.is_err());
    assert_eq!(error, Err(SmilesError::UnexpectedLeftBracket));
}

#[test]
fn test_unclosed_brackets() {
    let error = TokenIter::from("[Co+3").collect::<Result<Vec<_>, SmilesError>>();
    assert!(error.is_err());
    assert_eq!(error, Err(SmilesError::UnclosedBracket));
}

#[test]
fn test_smiles_tokens_water() {
    // "[OH2]" now produces ONE Token::BracketedAtom(...)
    let water_line = SMILES_WITH_BRACKETS[0];
    let tokens = TokenIter::from(water_line)
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_else(|_| panic!("Failed to parse {water_line}"));

    assert_eq!(tokens.len(), 1);

    match &tokens[0] {
        Token::BracketedAtom(b) => {
            // validate the important semantics
            assert_eq!(b.element(), Some(elements_rs::Element::O));
            assert!(!b.aromatic());
            assert_eq!(b.hydrogen_count(), Some(2));
            assert_eq!(b.charge_value(), 0);
            assert_eq!(b.class(), 0);
            assert_eq!(b.chiral(), None);
            assert_eq!(b.isotope_mass_number(), None);
        }
        other => panic!("Expected Token::BracketedAtom, got {other:?}"),
    }
}
