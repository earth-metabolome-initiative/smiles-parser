#![no_main]

use libfuzzer_sys::fuzz_target;
use smiles_parser::bond::Bond;
use smiles_parser::prelude::{AromaticityPolicy, AromaticityStatus, Smiles};

const MAX_FUZZ_SMILES_BYTES: usize = 512;
const MAX_FUZZ_ATOMS: usize = 512;

fn has_source_aromatic_labels(smiles: &Smiles) -> bool {
    smiles.nodes().iter().any(|atom| atom.aromatic())
        || smiles
            .nodes()
            .iter()
            .enumerate()
            .any(|(atom_id, _)| smiles.edges_for_node(atom_id).iter().any(|edge| edge.bond() == Bond::Aromatic))
}

fn assert_policy_roundtrip(smiles: &Smiles, policy: AromaticityPolicy) {
    let perception = match smiles.perceive_aromaticity_for(policy) {
        Ok(perception) => perception,
        Err(_) => return,
    };

    let preserved = perception
        .kekulize()
        .expect("preserve-source kekulization should succeed after perception");
    assert_eq!(
        &preserved,
        smiles,
        "preserve-source kekulization should restore the original graph for policy {policy:?}"
    );

    if perception.status() == AromaticityStatus::Unsupported {
        return;
    }

    let expected_atom_ids = perception.assignment().atom_ids().to_vec();
    let expected_bond_edges = perception.assignment().bond_edges().to_vec();

    let Ok(standalone) = perception.kekulize_standalone() else {
        return;
    };
    let reperceived = standalone
        .perceive_aromaticity_for(policy)
        .expect("reperception should succeed after standalone kekulization");

    assert_eq!(
        reperceived.assignment().atom_ids(),
        expected_atom_ids,
        "standalone kekulization changed aromatic atom ids for policy {policy:?}"
    );
    assert_eq!(
        reperceived.assignment().bond_edges(),
        expected_bond_edges,
        "standalone kekulization changed aromatic bond edges for policy {policy:?}"
    );
}

fuzz_target!(|data: &[u8]| {
    if data.len() > MAX_FUZZ_SMILES_BYTES {
        return;
    }

    let Ok(data) = core::str::from_utf8(data) else {
        return;
    };
    let Ok(smiles) = data.parse::<Smiles>() else {
        return;
    };
    if smiles.nodes().len() > MAX_FUZZ_ATOMS {
        return;
    }
    if has_source_aromatic_labels(&smiles) {
        return;
    }

    for policy in [
        AromaticityPolicy::RdkitDefault,
        AromaticityPolicy::RdkitSimple,
        AromaticityPolicy::RdkitMdl,
    ] {
        assert_policy_roundtrip(&smiles, policy);
    }
});
