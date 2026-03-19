//! Rendering-focused tests for SMILES output.

use std::str::FromStr;

use smiles_parser::{bond::Bond, smiles::Smiles};

const RENDER_EXACT: &[&str] = &[
    "C1=CC=CC=C1",
    "c1ccccc1",
    "c1ccccc1*",
    "CCO",
    "C#N",
    "[Na+].[Cl-]",
    "FC(Br)(Cl)F",
    "C1CCCC=1",
    "C1:C:C:C:C:C1",
    "c1ccccc1c2ccccc2",
    r"c1=c\c2ccc(cc2)-c2ccccc2-c2ccc/1cc2",
];

const RENDER_ROUND_TRIP: &[&str] = &[
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
    "CCN1C[C@]2(COC)CC[C@H](O)[C@@]34[C@@H]5C[C@H]6[C@H](OC)[C@@H]5[C@](O)(C[C@@H]6OC)[C@@](O)([C@@H](OC)[C@H]23)[C@@H]14",
    "C[C@@H]1C[C@@]2(O[C@H]2C)C(=O)O[C@@H]2CCN(C)C/C=C(/COC(=O)[C@]1(C)O)C2=O",
    "CC=C(C)C1=C(Cl)C(O)=C(C)C2=C1OC1=CC(O)=C(Cl)C(C)=C1C(=O)O2",
    "CC1=C[C@H](O)CC(C)(C)[C@H]1/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(\\C)CO",
    "c1ccccc1*",
    "COC1(C23C14C5=C6C7=C8C5=C9C1=C5C%10=C%11C%12=C%13C%10=C%10C1=C8C1=C%10C8=C%10C%14=C%15C%16=C%17C(=C%12C%12=C%17C%17=C%18C%16=C%16C%15=C%15C%10=C1C7=C%15C1=C%16C(=C%18C7=C2C2=C%10C(=C5C9=C42)C%11=C%12C%10=C%177)C3=C16)C%14=C%138)OC",
    "C1CCCC=1",
];

fn parse_or_panic(s: &str) -> Smiles {
    Smiles::from_str(s).unwrap_or_else(|e| panic!("Failed to parse:\n{}", e.render(s)))
}

fn bond_count(smiles: &Smiles, bond: Bond) -> usize {
    smiles.edges().iter().filter(|e| e.bond() == bond).count()
}

#[test]
fn test_render_exact_curated_cases() {
    for &original in RENDER_EXACT {
        let parsed = parse_or_panic(original);
        let rendered = parsed.to_string();

        assert_eq!(
            rendered, original,
            "Rendered text changed.\nOriginal:\n{original}\nRendered:\n{rendered}"
        );
    }
}

#[test]
fn test_render_round_trip_preserves_graph_shape_and_bond_counts() {
    for &original in RENDER_ROUND_TRIP {
        let parsed = parse_or_panic(original);
        let rendered = parsed.to_string();
        let reparsed = Smiles::from_str(&rendered).unwrap_or_else(|e| {
            panic!(
                "Failed to parse rendered SMILES.\nOriginal:\n{original}\nRendered:\n{rendered}\n{}",
                e.render(&rendered)
            )
        });

        assert_eq!(
            parsed.nodes().len(),
            reparsed.nodes().len(),
            "Node count changed after render round trip.\nOriginal:\n{original}\nRendered:\n{rendered}"
        );

        assert_eq!(
            parsed.edges().len(),
            reparsed.edges().len(),
            "Edge count changed after render round trip.\nOriginal:\n{original}\nRendered:\n{rendered}"
        );

        let parsed_aromatic = parsed.nodes().iter().filter(|n| n.atom().aromatic()).count();
        let reparsed_aromatic = reparsed.nodes().iter().filter(|n| n.atom().aromatic()).count();

        assert_eq!(
            parsed_aromatic, reparsed_aromatic,
            "Aromatic atom count changed after render round trip.\nOriginal:\n{original}\nRendered:\n{rendered}"
        );

        for bond in [
            Bond::Single,
            Bond::Double,
            Bond::Triple,
            Bond::Quadruple,
            Bond::Aromatic,
            Bond::Up,
            Bond::Down,
        ] {
            assert_eq!(
                bond_count(&parsed, bond),
                bond_count(&reparsed, bond),
                "Bond count for {bond:?} changed after render round trip.\nOriginal:\n{original}\nRendered:\n{rendered}"
            );
        }
    }
}

#[test]
fn test_render_preserves_explicit_closure_double_bond_text() {
    let original = "C1CCCC=1";
    let parsed = parse_or_panic(original);
    let rendered = parsed.to_string();

    assert_eq!(rendered, original);
}

#[test]
fn test_render_preserves_disconnected_components_text() {
    let original = "[Na+].[Cl-]";
    let parsed = parse_or_panic(original);
    let rendered = parsed.to_string();

    assert_eq!(rendered, original);
    assert_eq!(rendered.matches('.').count(), 1);
}

#[test]
fn test_render_preserves_branching_text() {
    let original = "FC(Br)(Cl)F";
    let parsed = parse_or_panic(original);
    let rendered = parsed.to_string();

    assert_eq!(rendered, original);
    assert_eq!(rendered.matches('(').count(), 2);
    assert_eq!(rendered.matches(')').count(), 2);
}

#[test]
fn test_render_preserves_aromatic_ring_text() {
    let original = "c1ccccc1";
    let parsed = parse_or_panic(original);
    let rendered = parsed.to_string();

    assert_eq!(rendered, original);
}

#[test]
fn test_render_preserves_ring_digits_on_benzene() {
    let original = "C1=CC=CC=C1";
    let parsed = parse_or_panic(original);
    let rendered = parsed.to_string();

    assert_eq!(rendered, original);
    assert_eq!(rendered.matches('1').count(), 2);
}
