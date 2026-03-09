//! Tests of the parser module for several corner cases.

use smiles_parser::{errors::SmilesErrorWithSpan, parser::token_iter::TokenIter};
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
    "COC1(C23C14C5=C6C7=C8C5=C9C1=C5C%10=C%11C%12=C%13C%10=C%10C1=C8C1=C%10C8=C%10C%14=C%15C%16=C%17C(=C%12C%12=C%17C%17=C%18C%16=C%16C%15=C%15C%10=C1C7=C%15C1=C%16C(=C%18C7=C2C2=C%10C(=C5C9=C42)C%11=C%12C%10=C%177)C3=C16)C%14=C%138)OC",
];

#[test]
fn test_tokenizer() {
    for &s in SMILES_STR {
        let _tokens = TokenIter::from(s)
            .collect::<Result<Vec<_>, SmilesErrorWithSpan>>()
            .unwrap_or_else(|e| panic!("Failed to tokenize:\n{}", e.render(s)));
    }
}

#[test]
fn test_parser_all_inputs() {
    use std::str::FromStr;

    use smiles_parser::smiles::Smiles;

    for &s in SMILES_STR {
        Smiles::from_str(s).unwrap_or_else(|e| panic!("Failed to parse:\n{}", e.render(s)));
    }
}

#[test]
fn test_parse_then_render_then_parse_all_inputs() {
    use std::str::FromStr;

    use smiles_parser::smiles::Smiles;

    for &original in SMILES_STR {
        let parsed = Smiles::from_str(original)
            .unwrap_or_else(|e| panic!("Failed to parse original:\n{}", e.render(original)));

        let rendered = parsed.to_string();

        let reparsed = Smiles::from_str(&rendered).unwrap_or_else(|e| {
            panic!(
                "Failed to parse rendered SMILES.\nOriginal: {original}\nRendered: {rendered}\n{}",
                e.render(&rendered)
            )
        });

        assert_eq!(
            parsed.nodes().len(),
            reparsed.nodes().len(),
            "Node count changed after round trip.\nOriginal: {original}\nRendered: {rendered}"
        );

        assert_eq!(
            parsed.edges().len(),
            reparsed.edges().len(),
            "Edge count changed after round trip.\nOriginal: {original}\nRendered: {rendered}"
        );

        let parsed_aromatic = parsed.nodes().iter().filter(|n| n.atom().aromatic()).count();
        let reparsed_aromatic = reparsed.nodes().iter().filter(|n| n.atom().aromatic()).count();

        assert_eq!(
            parsed_aromatic,
            reparsed_aromatic,
            "Aromatic atom count changed after round trip.\nOriginal: {original}\nRendered: {rendered}"
        );
    }
}

#[test]
fn test_parse_benzene_graph_shape() {
    // C1=CC=CC=C1
    use std::str::FromStr;

    use smiles_parser::smiles::Smiles;

    let line = SMILES_STR[0];
    let smiles =
        Smiles::from_str(line).unwrap_or_else(|e| panic!("Failed to parse:\n{}", e.render(line)));

    assert_eq!(smiles.nodes().len(), 6);
    assert_eq!(smiles.edges().len(), 6);
    assert!(smiles.nodes().iter().all(|n| !n.atom().aromatic()));
}

#[test]
fn test_parse_benzene_with_wildcard_graph_shape() {
    // c1ccccc1*
    use std::str::FromStr;

    use smiles_parser::smiles::Smiles;

    let line = SMILES_STR[31];
    let smiles =
        Smiles::from_str(line).unwrap_or_else(|e| panic!("Failed to parse:\n{}", e.render(line)));

    assert_eq!(smiles.nodes().len(), 7);
    assert_eq!(smiles.edges().len(), 7);
    let aromatic_count = smiles.nodes().iter().filter(|n| n.atom().aromatic()).count();
    assert_eq!(aromatic_count, 6);
}