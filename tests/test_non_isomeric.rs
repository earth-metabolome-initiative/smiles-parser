//! Non-isomeric normalization tests.
use smiles_parser::prelude::Smiles;

fn ni(s: &str) -> String {
    s.parse::<Smiles>().unwrap().non_isomeric().render()
}

#[test]
fn drops_isotopes() {
    assert!(!ni("[13CH3]O").contains("13"));
    assert_eq!(ni("[13CH3]O"), ni("[12CH3]O"));
}

#[test]
fn drops_tetrahedral_chirality() {
    let r = ni("N[C@@H](C)C(=O)O");
    assert!(!r.contains('@'), "got {r}");
    assert_eq!(ni("N[C@@H](C)C(=O)O"), ni("N[C@H](C)C(=O)O"));
}

#[test]
fn drops_double_bond_stereo() {
    let r = ni("F/C=C/F");
    assert!(!r.contains('/') && !r.contains('\\'), "got {r}");
    assert_eq!(ni("F/C=C/F"), ni("F/C=C\\F"));
}

#[test]
fn preserves_skeleton_and_is_idempotent() {
    for s in ["[13CH3]O", "N[C@@H](C)C(=O)O", "F/C=C/F", "c1ccccc1", "CCO"] {
        let once = s.parse::<Smiles>().unwrap().non_isomeric();
        assert_eq!(once.render(), once.non_isomeric().render(), "{s} not idempotent");
        // same heavy-atom count
        assert_eq!(once.nodes().len(), s.parse::<Smiles>().unwrap().nodes().len(), "{s}");
    }
}
