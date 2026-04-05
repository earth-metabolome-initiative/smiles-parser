//! Test for tokenizing square brackets

use std::str::FromStr;

use smiles_parser::{errors::SmilesError, smiles::Smiles};

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
    "CCN1C[C@]2(COC)CC[C@H](O)[C@@]34[C@@H]5C[C@H]6[C@H](OC)[C@@H]5[C@](O)(C[C@@H]6OC)[C@@](O)([C@@H](OC)[C@H]23)[C@@H]14",
    "C[C@@H]1C[C@@]2(O[C@H]2C)C(=O)O[C@@H]2CCN(C)C/C=C(/COC(=O)[C@]1(C)O)C2=O",
    "CC1=C[C@H](O)CC(C)(C)[C@H]1/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(\\C)CO",
];
#[test]
fn test_valid_square_brackets() {
    for &s in SMILES_WITH_BRACKETS {
        let _smiles = Smiles::from_str(s).unwrap_or_else(|_| panic!("Failed to parse {s}"));
    }
}

#[test]
fn test_unexpected_right_brackets() {
    // stray ']' becomes UnexpectedCharacter(']') with a span
    let err = Smiles::from_str("[Co+3]]").unwrap_err();

    assert_eq!(err.smiles_error(), SmilesError::UnexpectedCharacter(']'));
    assert_eq!(err.start(), 6);
    assert_eq!(err.end(), 7);
}

#[test]
fn test_unexpected_left_brackets() {
    let err = Smiles::from_str("[[Co+3]").unwrap_err();

    assert_eq!(err.smiles_error(), SmilesError::MissingElement);
    assert_eq!(err.start(), 0);
    assert_eq!(err.end(), 2);
}

#[test]
fn test_unclosed_brackets() {
    let err = Smiles::from_str("[Co+3").unwrap_err();

    assert_eq!(err.smiles_error(), SmilesError::UnclosedBracket);
    assert_eq!(err.start(), 0);
    assert_eq!(err.end(), 5);
}

#[test]
fn test_smiles_tokens_water() {
    let water_line = SMILES_WITH_BRACKETS[0];
    let smiles =
        Smiles::from_str(water_line).unwrap_or_else(|_| panic!("Failed to parse {water_line}"));

    assert_eq!(smiles.nodes().len(), 1);

    let atom = &smiles.nodes()[0];
    assert_eq!(atom.element(), Some(elements_rs::Element::O));
    assert!(!atom.aromatic());
    assert_eq!(atom.hydrogen_count(), 2);
    assert_eq!(atom.charge_value(), 0);
    assert_eq!(atom.class(), 0);
    assert_eq!(atom.chirality(), None);
    assert_eq!(atom.isotope_mass_number(), None);
}
