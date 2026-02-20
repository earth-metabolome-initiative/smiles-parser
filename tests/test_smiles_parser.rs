//! Tests of the parser module for several corner cases.

use elements_rs::Element;
use smiles_parser::{
    atom_symbol::AtomSymbol,
    errors::{SmilesError, SmilesErrorWithSpan},
    parser::token_iter::TokenIter,
    ring_num::RingNum,
    token::{Token, TokenWithSpan},
    unbracketed::UnbracketedAtom,
};
const SMILES_STR: &[&str] = &[
    "C1=CC=CC=C1",
    "[OH2]",
    "[Ti+4]",
    "[Co+3]",
    "CCO",
    "C#N",
    "[Ga+]$[As-]",
    "[Na+].[Cl-]",
    "C1CCCC2C1CCCC2",
    "C0CCCCC0C0CCCCC0",
    "C1:C:C:C:C:C1",
    "c1ccccc1",
    "c1ccccc1c2ccccc2",
    "COc(c1)cccc1C#N",
    "FC(Br)(Cl)F",
    "CC1CCC/C(C)=C1/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)/CCCC2(C)C",
    "N[C@@H](C)C(=O)O",
    "N[CH](C)C(=O)O",
    "NC(C)C(=O)O",
    "[14cH]1ccccc1",
    "[14c@H]1ccccc1",
    "[2H]C(Cl)(Cl)Cl",
    "[C@@H](C)(N)C(=O)O",
    "C[C@H](N)C(=O)O",
    "OC(=O)[C@@H](N)C",
    "[K+].C=C.Cl[Pt-](Cl)Cl.O",
    "[Ti+4]",
    "CCN1C[C@]2(COC)CC[C@H](O)[C@@]34[C@@H]5C[C@H]6[C@H](OC)[C@@H]5[C@](O)(C[C@@H]6OC)[C@@](O)([C@@H](OC)[C@H]23)[C@@H]14",
    "C[C@@H]1C[C@@]2(O[C@H]2C)C(=O)O[C@@H]2CCN(C)C/C=C(/COC(=O)[C@]1(C)O)C2=O",
    "CC=C(C)C1=C(Cl)C(O)=C(C)C2=C1OC1=CC(O)=C(Cl)C(C)=C1C(=O)O2",
    "CC1=C[C@H](O)CC(C)(C)[C@H]1/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(\\C)CO",
    "c1ccccc1*",
];

#[test]
fn test_tokenizer() -> Result<(), SmilesErrorWithSpan> {
    for &s in SMILES_STR {
        let _tokens = TokenIter::from(s)
            .collect::<Result<Vec<_>, SmilesErrorWithSpan>>()
            .unwrap_or_else(|e| panic!("Failed to tokenize:\n{}", e.render(s)));
    }
    Ok(())
}

#[test]
fn test_smiles_tokens_benzene() -> Result<(), SmilesError> {
    // C1=CC=CC=C1
    // (not aromatic because uppercase C)
    let c = UnbracketedAtom::new(AtomSymbol::Element(Element::C), false);

    let expected = vec![
        TokenWithSpan::new(Token::UnbracketedAtom(c), 0, 1),
        TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1)?), 1, 2),
        TokenWithSpan::new(Token::Bond(smiles_parser::bond::Bond::Double), 2, 3),
        TokenWithSpan::new(Token::UnbracketedAtom(c), 3, 4),
        TokenWithSpan::new(Token::UnbracketedAtom(c), 4, 5),
        TokenWithSpan::new(Token::Bond(smiles_parser::bond::Bond::Double), 5, 6),
        TokenWithSpan::new(Token::UnbracketedAtom(c), 6, 7),
        TokenWithSpan::new(Token::UnbracketedAtom(c), 7, 8),
        TokenWithSpan::new(Token::Bond(smiles_parser::bond::Bond::Double), 8, 9),
        TokenWithSpan::new(Token::UnbracketedAtom(c), 9, 10),
        TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1)?), 10, 11),
    ];
    let line = SMILES_STR[0];
    let got = TokenIter::from(line)
        .collect::<Result<Vec<_>, SmilesErrorWithSpan>>()
        .unwrap_or_else(|e| panic!("Failed to tokenize:\n{}", e.render(line)));

    assert_eq!(expected, got);
    Ok(())
}

#[test]
fn test_smiles_tokens_benzene_with_wildcard() -> Result<(), SmilesError> {
    // c1ccccc1*
    let aromatic_c = UnbracketedAtom::new(AtomSymbol::Element(Element::C), true);
    let star = UnbracketedAtom::new(AtomSymbol::WildCard, false);

    let expected = vec![
        TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 0, 1),
        TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1)?), 1, 2),
        TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 2, 3),
        TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 3, 4),
        TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 4, 5),
        TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 5, 6),
        TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 6, 7),
        TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1)?), 7, 8),
        TokenWithSpan::new(Token::UnbracketedAtom(star), 8, 9),
    ];

    let line = SMILES_STR[31];
    let got = TokenIter::from(line)
        .collect::<Result<Vec<_>, SmilesErrorWithSpan>>()
        .unwrap_or_else(|e| panic!("Failed to tokenize:\n{}", e.render(line)));

    assert_eq!(expected, got);
    Ok(())
}
