//! Implicit hydrogen counting for parsed SMILES graphs.
//!
//! This module intentionally implements the local, syntax-level SMILES rules
//! for the graph exactly as it was parsed. It does not perform aromaticity
//! perception, Kekule normalization, valence cleanup, or any other sanitization
//! step before counting hydrogens.
//!
//! In practice that means:
//! - bracket atoms always contribute zero implicit hydrogens unless a future
//!   API explicitly chooses to model bracket defaults differently. This matches
//!   `OpenSMILES` bracket-`H0` semantics and raw `RDKit` behavior.
//! - unbracketed organic-subset atoms use normal-valence completion. This also
//!   matches raw `RDKit` for unsanitized SMILES input.
//! - lowercase aromatic atoms use a small SMILES-specific policy table. The
//!   current table is chosen to match raw `RDKit` for atoms that were already
//!   parsed as aromatic tokens.
//! - aromatic bonds are counted as bond order 1 because this is a raw parsed
//!   graph, not a normalized Lewis structure. This also matches raw `RDKit`'s
//!   property-cache behavior before aromaticity perception or Kekule
//!   normalization.
//!
//! The current behavior is designed to match raw `RDKit` property-cache
//! semantics (`MolFromSmiles(..., sanitize=False)` followed by
//! `UpdatePropertyCache` with `strict=False`) for SMILES-as-written input.

use alloc::vec::Vec;

use elements_rs::{AllowedValences, ChargedValences, Element};
use geometric_traits::traits::SparseValuedMatrix2DRef;

use super::Smiles;
use crate::{
    atom::{Atom, AtomSyntax, atom_symbol::AtomSymbol},
    bond::Bond,
};

impl Smiles {
    /// Returns the per-atom implicit hydrogen counts in node order.
    ///
    /// This method is deliberately narrow: it answers the question "how many
    /// implicit hydrogens does this parsed SMILES tokenization imply right
    /// now?"
    ///
    /// It does not try to reinterpret the graph chemically. In particular, it
    /// does not:
    /// - perceive aromaticity for Kekule inputs
    /// - reassign bond orders
    /// - normalize charges or hypervalent patterns
    /// - validate that a lowercase aromatic spelling is chemically realizable
    ///
    /// That keeps this API cheap, deterministic, and aligned with SMILES-as-
    /// written semantics.
    ///
    /// Wherever this crate had to choose among multiple plausible local
    /// behaviors, the current implementation prefers the behavior observed from
    /// raw `RDKit` (`sanitize=False`, then
    /// `UpdatePropertyCache(strict=False)`).
    #[inline]
    #[must_use]
    pub fn implicit_hydrogen_counts(&self) -> Vec<u8> {
        if let Some(cached) = &self.implicit_hydrogen_cache {
            return cached.clone();
        }
        self.nodes()
            .iter()
            .enumerate()
            .map(|(node_id, node)| implicit_hydrogens_for_node(self, node_id, node))
            .collect()
    }

    /// Returns the implicit hydrogen count for a node id, if present.
    ///
    /// This uses the same raw parsed-graph semantics as
    /// [`Smiles::implicit_hydrogen_counts`].
    #[inline]
    #[must_use]
    pub fn implicit_hydrogen_count(&self, id: usize) -> Option<u8> {
        if let Some(cached) = &self.implicit_hydrogen_cache {
            return cached.get(id).copied();
        }
        self.node_by_id(id).map(|node| implicit_hydrogens_for_node(self, id, node))
    }
}

/// Computes the implicit hydrogen count for a single node using only local
/// graph information.
///
/// For bracket atoms, returning `0` is both the SMILES rule and the behavior
/// observed from raw `RDKit`.
#[inline]
fn implicit_hydrogens_for_node(smiles: &Smiles, node_id: usize, node: &Atom) -> u8 {
    let explicit_valence = explicit_valence(smiles, node_id);
    match node.syntax() {
        AtomSyntax::Bracket => 0,
        AtomSyntax::OrganicSubset => {
            match node.symbol() {
                AtomSymbol::WildCard => 0,
                AtomSymbol::Element(element) => {
                    if node.aromatic() {
                        aromatic_implicit_hydrogens(element, explicit_valence)
                    } else {
                        aliphatic_implicit_hydrogens(element, explicit_valence)
                    }
                }
            }
        }
    }
}

/// Returns the raw explicit valence contribution from the parsed graph.
///
/// Only bond orders are included here. Explicit bracket hydrogens are not added
/// because this helper is only used for unbracketed atoms; bracket atoms return
/// zero implicit hydrogens before valence completion is considered.
///
/// This split also mirrors raw `RDKit`: bracket hydrogens stay explicit instead
/// of being folded back into a later implicit-hydrogen completion step.
#[inline]
fn explicit_valence(smiles: &Smiles, node_id: usize) -> u8 {
    smiles.bond_matrix().sparse_row_values_ref(node_id).map(|entry| bond_order(entry.bond())).sum()
}

/// Maps parsed bond syntax to the local bond-order contribution used for raw
/// implicit-hydrogen counting.
///
/// Aromatic bonds contribute `1` here because this module works on the parsed
/// aromatic graph directly instead of first assigning a Kekule form.
/// That choice matches raw `RDKit` on unsanitized molecules before any aromatic
/// normalization pass is applied.
#[inline]
fn bond_order(bond: Bond) -> u8 {
    match bond {
        Bond::Single | Bond::Up | Bond::Down | Bond::Aromatic => 1,
        Bond::Double => 2,
        Bond::Triple => 3,
        Bond::Quadruple => 4,
    }
}

/// Applies normal-valence completion for an unbracketed aliphatic atom.
#[inline]
fn aliphatic_implicit_hydrogens(element: Element, explicit_valence: u8) -> u8 {
    target_valence(element, explicit_valence)
        .map_or(0, |target| target.saturating_sub(explicit_valence))
}

/// Selects the first compatible target valence for an unbracketed atom.
///
/// Most elements delegate to `elements-rs`, but a few halogens are handled
/// explicitly in order to stay aligned with raw `RDKit` behavior on unsanitized
/// SMILES input:
/// - `F`, `Cl`, and `Br` cap at valence 1 for implicit-hydrogen purposes
/// - neutral `I` advances through the sequence 1, 3, 5
///
/// These overrides are intentionally documented here because they are not just
/// periodic-table defaults; they are compatibility choices made to match raw
/// `RDKit`.
#[inline]
fn target_valence(element: Element, explicit_valence: u8) -> Option<u8> {
    if matches!(element, Element::F | Element::Cl | Element::Br) {
        return (explicit_valence <= 1).then_some(1);
    }
    if element == Element::I {
        return [1, 3, 5].into_iter().find(|candidate| *candidate >= explicit_valence);
    }

    let charged_valences = element.valences_at_charge(0);
    if charged_valences.is_empty() {
        let allowed_valences = element.allowed_valences();
        return allowed_valences.iter().copied().find(|candidate| *candidate >= explicit_valence);
    }

    charged_valences.iter().copied().find(|candidate| *candidate >= explicit_valence)
}

/// Applies the SMILES-specific aromatic defaults used by this crate.
///
/// This is intentionally a small policy table, not a general aromaticity model.
/// It only answers the local implicit-hydrogen question for atoms that were
/// already parsed as aromatic.
///
/// The table is also chosen to agree with raw `RDKit` for aromatic tokens such
/// as `c`, `n`, `o`, `s`, `b`, `p`, and bracketed aromatic atoms that default
/// to `H0` unless an explicit hydrogen count is written.
#[inline]
fn aromatic_implicit_hydrogens(element: Element, explicit_valence: u8) -> u8 {
    if element == Element::C { 3_u8.saturating_sub(explicit_valence) } else { 0 }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use elements_rs::{AllowedValences, ChargedValences, Element};

    use super::{
        Smiles, aromatic_implicit_hydrogens, bond_order, explicit_valence, target_valence,
    };
    use crate::bond::Bond;

    #[test]
    fn bracket_atoms_never_gain_implicit_hydrogens() {
        let smiles = Smiles::from_str("[CH3][OH-][NH4+]").unwrap();
        assert_eq!(smiles.implicit_hydrogen_counts(), vec![0, 0, 0]);
    }

    #[test]
    fn unbracketed_halogens_stop_at_valence_one() {
        let smiles = Smiles::from_str("FClF").unwrap();
        assert_eq!(smiles.implicit_hydrogen_counts(), vec![0, 0, 0]);
    }

    #[test]
    fn wildcard_atoms_never_gain_implicit_hydrogens() {
        let smiles = Smiles::from_str("*").unwrap();
        assert_eq!(smiles.implicit_hydrogen_counts(), vec![0]);
    }

    #[test]
    fn aromatic_defaults_match_current_policy() {
        let smiles = Smiles::from_str("c1ccccc1").unwrap();
        assert_eq!(smiles.implicit_hydrogen_counts(), vec![1, 1, 1, 1, 1, 1]);

        let pyridine = Smiles::from_str("n1ccccc1").unwrap();
        assert_eq!(pyridine.implicit_hydrogen_counts(), vec![0, 1, 1, 1, 1, 1]);
    }

    #[test]
    fn implicit_hydrogen_count_accesses_node_by_id() {
        let smiles = Smiles::from_str("CO").unwrap();
        assert_eq!(smiles.implicit_hydrogen_count(0), Some(3));
        assert_eq!(smiles.implicit_hydrogen_count(1), Some(1));
        assert_eq!(smiles.implicit_hydrogen_count(99), None);
    }

    #[test]
    fn explicit_valence_counts_quadruple_bonds() {
        let smiles = Smiles::from_str("C$C").unwrap();
        assert_eq!(explicit_valence(&smiles, 0), 4);
        assert_eq!(explicit_valence(&smiles, 1), 4);
        assert_eq!(smiles.implicit_hydrogen_counts(), vec![0, 0]);
    }

    #[test]
    fn bond_order_maps_double_and_triple_bonds() {
        assert_eq!(bond_order(Bond::Double), 2);
        assert_eq!(bond_order(Bond::Triple), 3);
    }

    #[test]
    fn target_valence_uses_neutral_iodine_progression() {
        assert_eq!(target_valence(Element::I, 2), Some(3));
        assert_eq!(target_valence(Element::I, 4), Some(5));
    }

    #[test]
    fn target_valence_falls_back_to_allowed_valences_when_neutral_table_is_empty() {
        assert_eq!(Element::Xe.valences_at_charge(0), &[] as &[u8]);
        assert_eq!(Element::Xe.allowed_valences(), &[0]);
        assert_eq!(target_valence(Element::Xe, 0), Some(0));
    }

    #[test]
    fn aromatic_fallback_defaults_unknown_cases_to_zero() {
        assert_eq!(aromatic_implicit_hydrogens(Element::Se, 0), 0);
    }
}
