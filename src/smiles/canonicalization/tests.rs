use alloc::vec::Vec;
use core::str::FromStr;

use geometric_traits::traits::SparseValuedMatrixRef;

use super::{
    Smiles, canonicalization_state_key, remap_parsed_stereo_neighbors_row,
    support::{
        assert_canonicalization_invariants, assert_distinct_canonicalization_state,
        hidden_bond_order_digest, permute_smiles, same_canonicalization_state,
    },
};
use crate::smiles::StereoNeighbor;

#[test]
fn canonical_labeling_of_empty_graph_is_empty() {
    let labeling = Smiles::new().canonical_labeling();
    assert!(labeling.order().is_empty());
    assert!(labeling.new_index_of_old_node().is_empty());
    assert!(Smiles::new().is_canonical());
}

#[test]
fn canonical_labeling_inverse_matches_order() {
    let smiles = Smiles::from_str("NCC(=O)O").unwrap();
    let labeling = smiles.canonical_labeling();

    assert_eq!(labeling.order().len(), smiles.nodes().len());
    let mut seen = labeling.order().to_vec();
    seen.sort_unstable();
    assert_eq!(seen, (0..smiles.nodes().len()).collect::<Vec<_>>());

    for (new_index, old_node) in labeling.order().iter().copied().enumerate() {
        assert_eq!(labeling.new_index_of_old_node()[old_node], new_index);
    }
}

#[test]
fn canonical_labeling_is_deterministic() {
    let smiles = Smiles::from_str("C1CC(C2CC2)CC1").unwrap();
    let first = smiles.canonical_labeling();
    let second = smiles.canonical_labeling();

    assert_eq!(first, second);
}

#[test]
fn canonicalize_is_idempotent_for_linear_case() {
    let smiles = Smiles::from_str("OCC").unwrap();
    let once = smiles.canonicalize();
    let twice = once.canonicalize();

    assert_eq!(once, twice);
    assert!(once.is_canonical());
}

#[test]
fn canonicalize_converges_equivalent_linear_spellings() {
    let left = Smiles::from_str("CCO").unwrap().canonicalize();
    let right = Smiles::from_str("OCC").unwrap().canonicalize();

    assert_eq!(left, right);
}

#[test]
fn canonicalize_converges_disconnected_component_permutations() {
    let left = Smiles::from_str("CC.O").unwrap().canonicalize();
    let right = Smiles::from_str("O.CC").unwrap().canonicalize();

    assert_eq!(left, right);
}

#[test]
fn canonicalize_converges_permuted_symmetric_cage_graph() {
    let original = Smiles::from_str("C12C3C4C1C5C2C3C45").unwrap();
    let permuted = permute_smiles(&original, &[7, 3, 0, 5, 2, 6, 1, 4]);

    assert_eq!(original.canonicalize(), permuted.canonicalize());
}

#[test]
fn canonicalize_converges_permuted_bridged_graph() {
    let original = Smiles::from_str("C1CC2CCC1C2").unwrap();
    let permuted = permute_smiles(&original, &[6, 2, 4, 0, 5, 1, 3]);

    assert_eq!(original.canonicalize(), permuted.canonicalize());
}

#[test]
fn canonicalize_converges_permuted_fused_graph() {
    let original = Smiles::from_str("c1cccc2ccccc12").unwrap();
    let permuted = permute_smiles(&original, &[9, 4, 0, 7, 2, 8, 1, 5, 3, 6]);

    assert_eq!(original.canonicalize(), permuted.canonicalize());
}

#[test]
fn canonicalize_marks_result_as_canonical() {
    let original = Smiles::from_str("OCC").unwrap();
    let canonicalized = original.canonicalize();

    assert!(!original.is_canonical());
    assert!(canonicalized.is_canonical());
}

#[test]
fn canonicalize_remaps_parsed_stereo_neighbors() {
    let smiles = Smiles::from_str("O[C@H](F)C").unwrap();
    let stereo_normalized = smiles.canonicalization_normal_form().stereo_normal_form();
    let labeling = stereo_normalized.exact_canonical_labeling();
    let canonicalized = smiles.canonicalize();

    let expected: Vec<Vec<StereoNeighbor>> = labeling
        .order()
        .iter()
        .copied()
        .map(|old_node| {
            remap_parsed_stereo_neighbors_row(
                &stereo_normalized,
                old_node,
                labeling.new_index_of_old_node(),
            )
        })
        .collect();

    let actual = (0..canonicalized.nodes().len())
        .map(|node_id| canonicalized.parsed_stereo_neighbors(node_id))
        .collect::<Vec<_>>();

    assert_eq!(actual, expected);
}

#[test]
fn canonicalize_remaps_implicit_hydrogen_cache_and_clears_provenance() {
    let aromaticized = Smiles::from_str("C1=CC=CC=C1")
        .unwrap()
        .perceive_aromaticity()
        .unwrap()
        .into_aromaticized();
    let labeling = aromaticized.canonical_labeling();
    let expected_cache = aromaticized
        .implicit_hydrogen_counts()
        .into_iter()
        .enumerate()
        .map(|(old_node, count)| (labeling.new_index_of_old_node()[old_node], count))
        .collect::<Vec<_>>();

    let canonicalized = aromaticized.canonicalize();
    let actual_cache =
        canonicalized.implicit_hydrogen_counts().into_iter().enumerate().collect::<Vec<_>>();

    let mut expected_cache = expected_cache;
    expected_cache.sort_unstable_by_key(|&(new_node, _)| new_node);

    assert_eq!(actual_cache, expected_cache);
    assert!(canonicalized.kekulization_source.is_none());
    assert!(
        canonicalized
            .bond_matrix()
            .sparse_entries()
            .all(|((row, column), entry)| row >= column || entry.ring_num().is_none())
    );
    assert_canonicalization_invariants(&aromaticized);
}

#[test]
fn canonicalize_preserves_isotope_distinction() {
    let natural = Smiles::from_str("CC").unwrap().canonicalize();
    let isotopic = Smiles::from_str("[13CH3]C").unwrap().canonicalize();

    assert_ne!(natural, isotopic);
}

#[test]
fn canonicalize_converges_benzene_aromatic_and_kekule_forms() {
    let aromatic = Smiles::from_str("c1ccccc1").unwrap().canonicalize();
    let kekule = Smiles::from_str("C1=CC=CC=C1").unwrap().canonicalize();

    assert_eq!(aromatic, kekule);
}

#[test]
fn canonicalize_converges_source_backed_aromatic_equivalence_groups() {
    let groups = [
        &["FC1=CC=CC=C1O", "FC1=C(O)C=CC=C1", "Fc1ccccc1O"][..],
        &["n1ccccc1", "N1=CC=CC=C1"][..],
        &["N1C=CC=C1", "[nH]1cccc1"][..],
        &["o1cccc1", "O1C=CC=C1"][..],
        &["s1cccc1", "S1C=CC=C1"][..],
        &["[cH-]1cccc1", "[CH-]1C=CC=C1"][..],
        &["c1cc2cccccc2c1", "C1=CC2=CC=CC=CC2=C1"][..],
        &["C1=CC=CN=C1", "c1cccnc1", "n1ccccc1"][..],
        &["[H]c2c([H])c(c1c(nc(n1([H]))C(F)(F)F)c2Cl)Cl", "Clc1ccc(Cl)c2[nH]c([nH0]c21)C(F)(F)F"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            assert_eq!(&canonicalized[0], other, "group did not converge: {group:?}");
        }
    }
}

#[test]
fn canonicalize_converges_additional_source_backed_connectivity_groups() {
    let groups = [
        &["OCC", "CCO", "C-C-O", "C(O)C"][..],
        &["C=CCC=CCO", "C=C-C-C=C-C-O", "OCC=CCC=C"][..],
        &["Oc1ccccc1", "c1ccccc1O", "c1(O)ccccc1", "c1(ccccc1)O"][..],
        &["C1=CCCCC1", "C=1CCCCC1", "C1CCCCC=1", "C=1CCCCC=1"][..],
        &["C1CCCCC1C1CCCCC1", "C1CCCCC1C2CCCCC2"][..],
        &["C12CCCCC1CCCC2", "C1CC2CCCCC2CC1"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            assert_eq!(&canonicalized[0], other, "group did not converge: {group:?}");
        }
    }
}

#[test]
fn canonicalize_converges_standard_form_bracket_and_explicit_hydrogen_groups() {
    let groups = [
        &["CCO", "[CH3]CO", "[CH3][CH2][OH]", "[H]OCC"][..],
        &["CC(=O)O", "[CH3]C(=O)O"][..],
        &["c1ccccc1", "[cH]1[cH][cH][cH][cH][cH]1"][..],
        &["Oc1ccccc1", "[OH][c]1[cH][cH][cH][cH][cH]1", "[H]Oc1ccccc1"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            same_canonicalization_state(&canonicalized[0], other);
        }
    }
}

#[test]
fn canonicalize_converges_additional_source_backed_real_world_groups() {
    let groups = [
        &["CC1=C(C=CO1)C(=O)NC2=CC=CC=C2", "Cc1occc1C(=O)Nc2ccccc2"][..],
        &[
            "CC1COC(Cn2cncn2)(c2ccc(Oc3ccc(Cl)cc3)cc2Cl)O1",
            "n1cnn(c1)CC2(OC(C)CO2)c3c(Cl)cc(cc3)Oc4ccc(Cl)cc4",
        ][..],
        &[
            "CC(C)C(Nc1ccc(C(F)(F)F)cc1Cl)C(=O)OC(C#N)c1cccc(Oc2ccccc2)c1",
            "FC(F)(F)c1cc(Cl)c(cc1)NC(C(C)(C))C(=O)OC(C#N)c2cccc(c2)Oc3ccccc3",
        ][..],
        &["NC(Cc1c[nH]c2ccccc12)C(=O)O", "OC(=O)C(N)C-c1c2ccccc2[nH]c1"][..],
        &["c1c(c23)ccc(c34)ccc4ccc2c1", "C1=CC2=C3C(=CC=C4C3=C1C=C4)C=C2"][..],
        &["S[Co@@](F)(Cl)(Br)(I)C=O", "O=C[Co@](F)(Cl)(Br)(I)S"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            same_canonicalization_state(&canonicalized[0], other);
        }
    }
}

#[test]
fn canonicalize_preserves_directional_bond_distinction() {
    let up = Smiles::from_str("F/C=C/F").unwrap().canonicalize();
    let down = Smiles::from_str("F/C=C\\F").unwrap().canonicalize();

    assert_ne!(up, down);
}

#[test]
fn canonicalize_converges_alkene_stereo_equivalence_groups() {
    let groups = [
        &["F/C=C/F", "F\\C=C\\F"][..],
        &["F/C=C\\F", "F\\C=C/F"][..],
        &["C/C=C/C", "C\\C=C\\C"][..],
        &["Cl/C=C\\Cl", "Cl\\C=C/Cl"][..],
        &["Cl/C=C/I", "Cl\\C=C\\I"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            same_canonicalization_state(&canonicalized[0], other);
        }
    }
}

#[test]
fn canonicalize_converges_atom_based_alkene_stereo_equivalence_groups() {
    let groups = [
        &["F[C@@H]=[C@H]F", "F[C@H]=[C@@H]F", "F/C=C/F", "[C@H](F)=[C@H]F"][..],
        &["F[C@H]=[C@H]F", "F[C@@H]=[C@@H]F", "F/C=C\\F", "[C@@H](F)=[C@H]F"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            same_canonicalization_state(&canonicalized[0], other);
        }
    }
}

#[test]
fn canonicalize_converges_tetrahedral_stereo_equivalence_groups() {
    let groups = [
        &[
            "N[C@](Br)(O)C",
            "Br[C@](O)(N)C",
            "O[C@](Br)(C)N",
            "Br[C@](C)(O)N",
            "C[C@](Br)(N)O",
            "Br[C@](N)(C)O",
            "C[C@@](Br)(O)N",
            "Br[C@@](N)(O)C",
            "[C@@](C)(Br)(O)N",
            "[C@@](Br)(N)(O)C",
        ][..],
        &["N[C@H](O)C", "O[C@H](C)N", "C[C@@H](O)N"][..],
        &["FC1C[C@](Br)(Cl)CCC1", "[C@]1(Br)(Cl)CCCC(F)C1"][..],
        &["C[C@H]1CCCCO1", "O1CCCC[C@@H]1C"][..],
        &["N[C@](Br)(O)C", "N[C@TH1](Br)(O)C"][..],
        &["N[C@@](Br)(O)C", "N[C@TH2](Br)(O)C"][..],
        &["[C@H](O)(N)C", "O[C@@H](N)C"][..],
        &["OC(Cl)=[C@]=C(C)F", "OC(Cl)=[C@AL1]=C(C)F"][..],
        &["OC(Cl)=[C@@]=C(C)F", "OC(Cl)=[C@AL2]=C(C)F"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            same_canonicalization_state(&canonicalized[0], other);
        }
    }
}

#[test]
fn canonicalize_clears_achiral_tetrahedral_markup() {
    let group = ["BrC(Br)C", "Br[C@H](Br)C", "Br[C@@H](Br)C"];
    let canonicalized = group
        .iter()
        .map(|source| Smiles::from_str(source).unwrap().canonicalize())
        .collect::<Vec<_>>();
    for other in &canonicalized[1..] {
        same_canonicalization_state(&canonicalized[0], other);
    }
}

#[test]
fn canonicalize_clears_non_stereogenic_alkene_markup() {
    let group = ["FC(F)=CF", "F/C(/F)=C/F", "F/C(/F)=C\\F"];
    let canonicalized = group
        .iter()
        .map(|source| Smiles::from_str(source).unwrap().canonicalize())
        .collect::<Vec<_>>();
    for other in &canonicalized[1..] {
        same_canonicalization_state(&canonicalized[0], other);
    }
}

#[test]
fn canonicalize_converges_non_tetrahedral_stereo_equivalence_groups() {
    let groups = [
        &["F[Po@SP1](Cl)(Br)I", "F[Po@SP2](Br)(Cl)I", "F[Po@SP3](Cl)(I)Br"][..],
        &["S[As@TB1](F)(Cl)(Br)N", "S[As@TB2](Br)(Cl)(F)N"][..],
        &["S[As@TB5](F)(N)(Cl)Br", "F[As@TB10](S)(Cl)(N)Br"][..],
        &["F[As@TB15](Cl)(S)(Br)N", "Br[As@TB20](Cl)(S)(F)N"][..],
        &["C[Co@OH1](F)(Cl)(Br)(I)S", "C[Co@OH8](F)(Br)(Cl)(I)S"][..],
        &["C[Co@](F)(Cl)(Br)(I)S", "F[Co@@](S)(I)(C)(Cl)Br"][..],
        &["S[Co@OH5](F)(I)(Cl)(C)Br", "Br[Co@OH9](C)(S)(Cl)(F)I"][..],
        &["Br[Co@OH12](Cl)(I)(F)(S)C", "Cl[Co@OH15](C)(Br)(F)(I)S"][..],
        &["Cl[Co@OH19](C)(I)(F)(S)Br", "I[Co@OH27](Cl)(Br)(F)(S)C"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            if let Err(_panic) = std::panic::catch_unwind(|| {
                same_canonicalization_state(&canonicalized[0], other);
            }) {
                panic!("group did not converge: {group:?}");
            }
        }
    }
}

#[test]
fn canonicalize_converges_source_backed_square_planar_groups() {
    let cis_group = [
        "[NH3][Pt@SP1]([NH3])(Cl)Cl",
        "[NH3][Pt@SP3]([NH3])(Cl)Cl",
        "[NH3][Pt@SP2](Cl)([NH3])Cl",
        "[NH3][Pt@SP1](Cl)(Cl)[NH3]",
    ];
    let trans_group =
        ["[NH3][Pt@SP2]([NH3])(Cl)Cl", "[NH3][Pt@SP1](Cl)([NH3])Cl", "[NH3][Pt@SP3](Cl)(Cl)[NH3]"];

    let cis_canonicalized = cis_group
        .iter()
        .map(|source| Smiles::from_str(source).unwrap().canonicalize())
        .collect::<Vec<_>>();
    for other in &cis_canonicalized[1..] {
        assert_eq!(&cis_canonicalized[0], other);
    }

    let trans_canonicalized = trans_group
        .iter()
        .map(|source| Smiles::from_str(source).unwrap().canonicalize())
        .collect::<Vec<_>>();
    for other in &trans_canonicalized[1..] {
        assert_eq!(&trans_canonicalized[0], other);
    }

    assert_distinct_canonicalization_state(&cis_canonicalized[0], &trans_canonicalized[0]);
}

#[test]
fn canonicalize_converges_square_planar_implicit_hydrogen_spellings() {
    let groups = [
        &["C[Pt@SP1H](Cl)F", "C[Pt@SP1]([H])(Cl)F"][..],
        &["C[Pt@SP1H2]Cl", "C[Pt@SP1]([H])([H])Cl"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            same_canonicalization_state(&canonicalized[0], other);
        }
    }
}

#[test]
fn canonicalize_converges_disconnected_ring_closure_stereo_variants() {
    let groups = [
        &["[C@@]1(Cl)(F)(I).Br1", "[C@@](Br)(Cl)(F)(I)"][..],
        &["[C@@](Cl)(F)(I)1.Br1", "[C@@](Cl)(F)(I)Br"][..],
    ];

    for group in groups {
        let canonicalized = group
            .iter()
            .map(|source| Smiles::from_str(source).unwrap().canonicalize())
            .collect::<Vec<_>>();
        for other in &canonicalized[1..] {
            same_canonicalization_state(&canonicalized[0], other);
        }
    }
}

#[test]
fn canonicalize_preserves_semantic_stereo_distinctions() {
    let alkene_e = Smiles::from_str("F/C=C/F").unwrap().canonicalize();
    let alkene_z = Smiles::from_str("F/C=C\\F").unwrap().canonicalize();
    let alkene_unspecified = Smiles::from_str("FC=CF").unwrap().canonicalize();
    let tetra_left = Smiles::from_str("N[C@](Br)(O)C").unwrap().canonicalize();
    let tetra_right = Smiles::from_str("N[C@@](Br)(O)C").unwrap().canonicalize();
    let allene_left = Smiles::from_str("OC(Cl)=[C@]=C(C)F").unwrap().canonicalize();
    let allene_right = Smiles::from_str("OC(Cl)=[C@@]=C(C)F").unwrap().canonicalize();

    assert_ne!(alkene_e, alkene_z);
    assert_ne!(alkene_e, alkene_unspecified);
    assert_ne!(tetra_left, tetra_right);
    assert_ne!(allene_left, allene_right);
}

#[test]
fn canonicalize_preserves_non_tetrahedral_stereo_distinctions() {
    let sp1 = Smiles::from_str("F[Po@SP1](Cl)(Br)I").unwrap().canonicalize();
    let sp2 = Smiles::from_str("F[Po@SP2](Cl)(Br)I").unwrap().canonicalize();
    let sp3 = Smiles::from_str("F[Po@SP3](Cl)(Br)I").unwrap().canonicalize();

    let tb1 = Smiles::from_str("S[As@TB1](F)(Cl)(Br)N").unwrap().canonicalize();
    let tb2 = Smiles::from_str("S[As@TB2](F)(Cl)(Br)N").unwrap().canonicalize();

    let oh1 = Smiles::from_str("C[Co@OH1](F)(Cl)(Br)(I)S").unwrap().canonicalize();
    let oh2 = Smiles::from_str("C[Co@OH2](F)(Cl)(Br)(I)S").unwrap().canonicalize();

    assert_distinct_canonicalization_state(&sp1, &sp2);
    assert_distinct_canonicalization_state(&sp1, &sp3);
    assert_distinct_canonicalization_state(&sp2, &sp3);
    assert_distinct_canonicalization_state(&tb1, &tb2);
    assert_distinct_canonicalization_state(&oh1, &oh2);
}

#[test]
fn canonicalize_handles_wildcard_quadruple_bond_regression() {
    let smiles = Smiles::from_str("*$*").unwrap();

    assert_canonicalization_invariants(&smiles);
}

#[test]
fn canonicalize_handles_aromatic_sulfur_chain_regression() {
    let smiles = Smiles::from_str("scCC").unwrap();

    assert_canonicalization_invariants(&smiles);
}

#[test]
fn canonicalize_handles_mixed_aromatic_phosphorus_wildcard_regression() {
    let smiles = Smiles::from_str("bsC3P1*sC3P1**p1C3P1**OC3p1C3Os*p1C3P1**OC3p1C3Os").unwrap();

    assert_canonicalization_invariants(&smiles);
}

#[test]
fn canonicalize_handles_polywildcard_aromatic_oxygen_regression() {
    let smiles =
        Smiles::from_str("OO*c8Oc8OOOO**c8OOO*c8OOOO*CO*c8OOOO*c8OOOO**c8OOO*c8OOOO*COOOO")
            .unwrap();

    assert_canonicalization_invariants(&smiles);
}

#[test]
fn canonicalize_reaches_fixed_point_for_pubchem_cid_10888631_regression() {
    let smiles =
        Smiles::from_str("C1C2[C@@H]3[C@H]4[C@H]5[C@H]6C5CC4[C@@H]6[C@H]2[C@@H]7[C@@H]3C71")
            .unwrap();

    let canonicalized = smiles.canonicalize();
    assert!(canonicalized.is_canonical(), "canonicalize() did not reach a canonical fixed point");

    let recanonicalized = canonicalized.canonicalize();
    same_canonicalization_state(&canonicalized, &recanonicalized);
}

#[test]
fn canonicalize_reaches_fixed_point_for_pubchem_cid_101510359_regression() {
    let smiles = Smiles::from_str("CC(C)(C)C1=CC(=CC(=C1)C2=C3C=CC(=C(C4=NC5=CC6=C(C=C(N6)C(=C7C=CC2=N7)C8=CC(=CC(=C8)C(C)(C)C)C(C)(C)C)C9=C1C=CC(=C(C2=NC(=C(C6=CC=C(N6)C(=C6C=CC9=N6)C6=CC(=CC(=C6)C(C)(C)C)C(C)(C)C)C6=CC(=CC(=C6)C(C)(C)C)C(C)(C)C)C=C2)C2=C6C=C7C(=CC(=N7)C(=C7C=CC(=C(C8=NC(=C(C(=C2)N6)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)C=C8)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)N7)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)C2=C6C=CC(=C(C5=C4)C4=NC(=C(C5=CC=C(N5)C(=C5C=CC2=N5)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)C=C4)N6)N1)C1=CC(=CC(=C1)C(C)(C)C)C(C)(C)C)N3)C(C)(C)C").unwrap();

    let canonicalized = smiles.canonicalize();
    let recanonicalized = canonicalized.canonicalize();

    same_canonicalization_state(&canonicalized, &recanonicalized);
}

#[test]
fn canonicalization_step_does_not_hide_bond_order_drift_for_pubchem_cid_101510359() {
    let smiles = Smiles::from_str("CC(C)(C)C1=CC(=CC(=C1)C2=C3C=CC(=C(C4=NC5=CC6=C(C=C(N6)C(=C7C=CC2=N7)C8=CC(=CC(=C8)C(C)(C)C)C(C)(C)C)C9=C1C=CC(=C(C2=NC(=C(C6=CC=C(N6)C(=C6C=CC9=N6)C6=CC(=CC(=C6)C(C)(C)C)C(C)(C)C)C6=CC(=CC(=C6)C(C)(C)C)C(C)(C)C)C=C2)C2=C6C=C7C(=CC(=N7)C(=C7C=CC(=C(C8=NC(=C(C(=C2)N6)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)C=C8)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)N7)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)C2=C6C=CC(=C(C5=C4)C4=NC(=C(C5=CC=C(N5)C(=C5C=CC2=N5)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)C2=CC(=CC(=C2)C(C)(C)C)C(C)(C)C)C=C4)N6)N1)C1=CC(=CC(=C1)C(C)(C)C)C(C)(C)C)N3)C(C)(C)C").unwrap();

    let entry = smiles.canonicalize();
    let direct_step = entry.canonicalization_step();

    assert_eq!(canonicalization_state_key(&entry), canonicalization_state_key(&direct_step));
    assert_eq!(hidden_bond_order_digest(&entry), hidden_bond_order_digest(&direct_step));
}

#[test]
fn canonicalize_converges_default_wildcard_parser_spellings() {
    let organic = Smiles::from_str("O*O").unwrap().canonicalize();
    let bracket = Smiles::from_str("O[*]O").unwrap().canonicalize();

    assert_eq!(organic, bracket);
}

#[test]
fn canonicalize_converges_permuted_graph_with_sidecars() {
    let original = Smiles::from_str("C1=CC=CC=C1[C@H](F)O")
        .unwrap()
        .perceive_aromaticity()
        .unwrap()
        .into_aromaticized();
    let permuted = permute_smiles(&original, &[8, 3, 6, 1, 7, 0, 4, 2, 5]);

    assert_eq!(original.canonicalize(), permuted.canonicalize());
}
