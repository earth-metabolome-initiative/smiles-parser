//! RDKit sanitization-reject regression cases.
//!
//! RDKit `2024.09.6` can parse these strings with `sanitize=False`, but rejects
//! them once sanitization runs because they are chemically invalid. This crate
//! currently accepts them because it does not run an equivalent validation
//! step.
//!
//! Expected work:
//! - add a post-parse validation pass, or route `FromStr` through one
//! - keep parser-grammar failures and chemical-validation failures distinct
//! - make the top-level API reject these inputs

use std::str::FromStr;

use smiles_parser::smiles::Smiles;

const RD_KIT_SANITIZATION_REJECTS: &[&str] = &[
    "[CH5]",
    "[OH3]",
    "[NH5+]",
    "C(F)(F)(F)(F)(F)F",
    "N(N)(N)(N)N",
    "C(C)(C)(C)(C)C",
    "O(O)(O)O",
    "[NH5]",
    "[OH4+]",
    "C(=O)(=O)=O",
    "N(=N)(=N)=N",
    "[BH5]",
    "[ClH2]",
];

#[test]
fn documented_rdkit_sanitization_rejects_should_fail_here_too() {
    let mut unexpected_accepts = Vec::new();

    for input in RD_KIT_SANITIZATION_REJECTS {
        match Smiles::from_str(&input) {
            Ok(smiles) => unexpected_accepts.push(format!("{input:?} parsed as {}", smiles)),
            Err(_) => {}
        }
    }

    assert!(
        unexpected_accepts.is_empty(),
        "documented RDKit sanitization rejects are still accepted here:\n{}",
        unexpected_accepts.join("\n")
    );
}
