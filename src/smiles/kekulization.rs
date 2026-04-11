//! Kekulization utilities for converting aromatic graphs to explicit
//! single/double-bond forms.

use alloc::vec::Vec;

use geometric_traits::traits::{
    Gabow1976, Matrix2D, SizedSparseMatrix2D, SizedSparseValuedMatrixRef, SparseValuedMatrix2DRef,
    SparseValuedMatrixRef,
};
use thiserror::Error;

use super::{BondEntry, BondMatrix, Smiles};
use crate::bond::Bond;

/// Error raised while converting aromatic bonds into a Kekule form.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum KekulizationError {
    /// No perfect matching exists for the aromatic candidate graph.
    #[error(
        "cannot assign a Kekule form for the aromatic candidate graph with {candidate_atom_count} candidate atoms"
    )]
    NoPerfectMatching {
        /// Number of candidate atoms that required localization.
        candidate_atom_count: usize,
    },
}

/// How kekulization should treat any preserved pre-aromatic source graph.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum KekulizationMode {
    /// Reuse the preserved pre-aromatic source graph when present.
    ///
    /// This is the natural inverse of
    /// [`Smiles::perceive_aromaticity`](super::Smiles::perceive_aromaticity)
    /// and keeps the original localized form when the aromatic graph was
    /// produced by this crate's aromaticity perception pipeline.
    PreserveSource,
    /// Ignore any preserved source graph and solve from the current aromatic
    /// graph alone.
    ///
    /// This is the mode to use when benchmarking or when the caller wants
    /// standalone kekulization semantics for an aromatic graph value.
    Standalone,
}

impl Smiles {
    /// Returns a localized Kekule form of the current graph.
    ///
    /// This is equivalent to calling
    /// [`Smiles::kekulize_with`](super::Smiles::kekulize_with) with
    /// [`KekulizationMode::PreserveSource`].
    ///
    /// # Examples
    ///
    /// ```rust
    /// use core::str::FromStr;
    ///
    /// use smiles_parser::prelude::Smiles;
    ///
    /// let original = Smiles::from_str("C1=CN=CN1").expect("valid Kekule imidazole");
    /// let perception = original.perceive_aromaticity().expect("perception should succeed");
    /// let restored = perception.aromaticized().kekulize().expect("kekulization should succeed");
    ///
    /// assert_eq!(restored, original);
    /// ```
    ///
    /// # Errors
    /// Returns [`KekulizationError::NoPerfectMatching`] if no valid localized
    /// assignment exists for the aromatic candidate graph.
    pub fn kekulize(&self) -> Result<Self, KekulizationError> {
        self.kekulize_with(KekulizationMode::PreserveSource)
    }

    /// Returns a localized Kekule form of the current graph while explicitly
    /// choosing how preserved pre-aromatic source graphs should be handled.
    ///
    /// Graphs with no aromatic bonds are returned unchanged.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use core::str::FromStr;
    ///
    /// use smiles_parser::prelude::{KekulizationMode, Smiles};
    ///
    /// let aromatic = Smiles::from_str("c1ccccc1").expect("valid aromatic benzene");
    /// let kekule = aromatic
    ///     .kekulize_with(KekulizationMode::Standalone)
    ///     .expect("standalone kekulization should succeed");
    /// let reperceived = kekule.perceive_aromaticity().expect("reperception should succeed");
    ///
    /// assert_eq!(reperceived.assignment().atom_ids(), &[0, 1, 2, 3, 4, 5]);
    /// assert_eq!(reperceived.assignment().bond_edges().len(), 6);
    /// ```
    ///
    /// # Errors
    /// Returns [`KekulizationError::NoPerfectMatching`] if no valid localized
    /// assignment exists for the aromatic candidate graph.
    pub fn kekulize_with(&self, mode: KekulizationMode) -> Result<Self, KekulizationError> {
        let aromatic_bonds = self
            .bond_matrix()
            .sparse_entries()
            .filter_map(|((row, column), entry)| {
                (row < column && entry.bond() == Bond::Aromatic).then_some(AromaticBond {
                    node_a: row,
                    node_b: column,
                    order: entry.order(),
                })
            })
            .collect::<Vec<_>>();

        if aromatic_bonds.is_empty() {
            return Ok(self.clone_without_kekulization_source());
        }

        if mode == KekulizationMode::PreserveSource
            && let Some(source) = &self.kekulization_source
        {
            return Ok((**source).clone());
        }

        let candidate_atom_ids = candidate_atom_ids(self);
        if has_unlocalizable_aromatic_component(
            self.nodes().len(),
            &aromatic_bonds,
            &candidate_atom_ids,
        ) {
            return Err(KekulizationError::NoPerfectMatching { candidate_atom_count: 0 });
        }

        let candidate_graph =
            KekulizationCandidateGraph::new(self, &aromatic_bonds, &candidate_atom_ids);
        if candidate_graph.graph.number_of_rows() == 0 {
            return Err(KekulizationError::NoPerfectMatching { candidate_atom_count: 0 });
        }
        let matched_orders = candidate_graph.matched_original_orders()?;
        let mut matched_order_flags = vec![false; self.number_of_bonds()];
        for matched_order in matched_orders {
            matched_order_flags[matched_order] = true;
        }

        let atom_nodes = self
            .atom_nodes
            .iter()
            .copied()
            .map(|atom| if atom.aromatic() { atom.with_aromatic(false) } else { atom })
            .collect::<Vec<_>>();

        let bond_matrix = BondMatrix::from_sorted_upper_triangular_entries(
            atom_nodes.len(),
            self.bond_matrix().sparse_entries().filter_map(|((row, column), entry)| {
                (row < column).then_some((
                    row,
                    column,
                    if entry.bond() == Bond::Aromatic {
                        if matched_order_flags[entry.order()] {
                            entry.with_bond(Bond::Double)
                        } else {
                            entry.with_bond(Bond::Single)
                        }
                    } else {
                        *entry
                    },
                ))
            }),
        )
        .unwrap_or_else(|_| unreachable!("existing bond matrix entries are already valid"));

        Ok(Self::from_bond_matrix_parts_with_implicit_hydrogen_cache(
            atom_nodes,
            bond_matrix,
            self.implicit_hydrogen_counts(),
        ))
    }

    /// Returns a localized Kekule form by solving from the current aromatic
    /// graph alone.
    ///
    /// This ignores any preserved pre-aromatic source graph.
    ///
    /// Aromatic graphs that do not retain enough local chemistry to be
    /// localized from their current state alone, such as some dummy-atom
    /// aromatic systems, return an error instead of guessing a bond pattern.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use core::str::FromStr;
    ///
    /// use smiles_parser::prelude::Smiles;
    ///
    /// let aromatic = Smiles::from_str("c1ccccc1").expect("valid aromatic benzene");
    /// let kekule = aromatic.kekulize_standalone().expect("standalone kekulization should succeed");
    ///
    /// assert!(kekule.to_string().contains('='));
    /// ```
    ///
    /// # Errors
    /// Returns [`KekulizationError::NoPerfectMatching`] if no valid localized
    /// assignment exists for the aromatic candidate graph.
    pub fn kekulize_standalone(&self) -> Result<Self, KekulizationError> {
        self.kekulize_with(KekulizationMode::Standalone)
    }
}

#[derive(Clone, Copy)]
struct AromaticBond {
    node_a: usize,
    node_b: usize,
    order: usize,
}

struct KekulizationCandidateGraph {
    graph: BondMatrix,
}

impl KekulizationCandidateGraph {
    fn new(smiles: &Smiles, aromatic_bonds: &[AromaticBond], local_atom_ids: &[usize]) -> Self {
        let mut original_to_local = vec![None; smiles.nodes().len()];
        for (local_id, atom_id) in local_atom_ids.iter().copied().enumerate() {
            original_to_local[atom_id] = Some(local_id);
        }

        let mut local_edges = Vec::with_capacity(aromatic_bonds.len());

        for aromatic_bond in aromatic_bonds {
            let AromaticBond { node_a, node_b, order } = *aromatic_bond;
            let (Some(local_a), Some(local_b)) =
                (original_to_local[node_a], original_to_local[node_b])
            else {
                continue;
            };
            local_edges.push((local_a, local_b, BondEntry::new(Bond::Single, None, order)));
        }

        if local_edges.len() > 1 {
            local_edges.sort_unstable_by_key(|(row, column, _)| (*row, *column));
        }

        let graph =
            BondMatrix::from_sorted_upper_triangular_entries(local_atom_ids.len(), local_edges)
                .unwrap_or_else(|_| {
                    unreachable!("candidate graph only adds unique non-self upper-triangular edges")
                });

        Self { graph }
    }

    fn matched_original_orders(&self) -> Result<Vec<usize>, KekulizationError> {
        let candidate_atom_count = self.graph.number_of_rows();
        if candidate_atom_count == 0 {
            return Ok(Vec::new());
        }
        if !candidate_atom_count.is_multiple_of(2) {
            return Err(KekulizationError::NoPerfectMatching { candidate_atom_count });
        }

        let matches = self.graph.gabow_1976();
        if matches.len() * 2 != candidate_atom_count {
            return Err(KekulizationError::NoPerfectMatching { candidate_atom_count });
        }

        Ok(matches
            .into_iter()
            .map(|(left, right)| {
                let rank =
                    self.graph.try_rank(left.min(right), left.max(right)).unwrap_or_else(|| {
                        unreachable!("matching edges always come from the candidate graph")
                    });
                self.graph.select_value_ref(rank).order()
            })
            .collect())
    }
}

fn has_unlocalizable_aromatic_component(
    atom_count: usize,
    aromatic_bonds: &[AromaticBond],
    candidate_atom_ids: &[usize],
) -> bool {
    let mut adjacency = vec![Vec::new(); atom_count];
    for bond in aromatic_bonds {
        adjacency[bond.node_a].push(bond.node_b);
        adjacency[bond.node_b].push(bond.node_a);
    }

    let mut candidate_flags = vec![false; atom_count];
    for atom_id in candidate_atom_ids {
        candidate_flags[*atom_id] = true;
    }

    let mut visited = vec![false; atom_count];
    let mut stack = Vec::new();

    for start in 0..atom_count {
        if visited[start] || adjacency[start].is_empty() {
            continue;
        }

        visited[start] = true;
        stack.push(start);
        let mut component_candidate_count = 0usize;

        while let Some(node) = stack.pop() {
            component_candidate_count += usize::from(candidate_flags[node]);
            for &neighbor in &adjacency[node] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }

        if component_candidate_count == 0 {
            return true;
        }
    }

    false
}

fn candidate_atom_ids(smiles: &Smiles) -> Vec<usize> {
    let implicit_hydrogen_counts = smiles.implicit_hydrogen_counts();
    smiles
        .nodes()
        .iter()
        .enumerate()
        .filter_map(|(atom_id, atom)| {
            needs_localized_double_bond(smiles, atom_id, atom, implicit_hydrogen_counts[atom_id])
                .then_some(atom_id)
        })
        .collect()
}

fn needs_localized_double_bond(
    smiles: &Smiles,
    atom_id: usize,
    atom: &crate::atom::Atom,
    implicit_hydrogen_count: u8,
) -> bool {
    if !atom.aromatic() {
        return false;
    }

    let Some(element) = atom.element() else {
        return false;
    };

    let mut preexisting_pi_bonds = 0usize;
    let mut degree = 0usize;
    for entry in smiles.bond_matrix().sparse_row_values_ref(atom_id) {
        degree += 1;
        match entry.bond().without_direction() {
            Bond::Double => preexisting_pi_bonds += 1,
            Bond::Triple | Bond::Quadruple => return false,
            _ => {}
        }
    }

    let total_hydrogens = usize::from(atom.hydrogen_count()) + usize::from(implicit_hydrogen_count);
    let valence = degree + total_hydrogens + preexisting_pi_bonds;

    match element {
        elements_rs::Element::B => {
            (atom.charge_value() == 0 && valence <= 2)
                || (atom.charge_value() == -1 && valence <= 3)
        }
        elements_rs::Element::C
        | elements_rs::Element::Si
        | elements_rs::Element::Ge
        | elements_rs::Element::Sn => {
            (atom.charge_value() == 0 && valence <= 3)
                || (atom.charge_value() == -1 && valence <= 2)
        }
        elements_rs::Element::N
        | elements_rs::Element::P
        | elements_rs::Element::As
        | elements_rs::Element::Sb => {
            (atom.charge_value() == 0 && (valence <= 2 || valence == 4))
                || (atom.charge_value() == 1 && valence <= 3)
        }
        elements_rs::Element::O
        | elements_rs::Element::S
        | elements_rs::Element::Se
        | elements_rs::Element::Te => {
            (atom.charge_value() == 0 && (valence <= 1 || valence == 3 || valence == 5))
                || (atom.charge_value() == 1 && (valence <= 2 || valence == 4))
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use alloc::{string::ToString, vec::Vec};
    use std::str::FromStr;

    use elements_rs::Element;
    use geometric_traits::traits::SparseValuedMatrixRef;

    use super::{
        AromaticBond, KekulizationCandidateGraph, KekulizationError, KekulizationMode, Smiles,
        candidate_atom_ids,
    };
    use crate::{
        atom::{Atom, atom_symbol::AtomSymbol},
        bond::{Bond, bond_edge::BondEdge},
        smiles::BondMatrixBuilder,
    };

    fn count_bonds(smiles: &Smiles, bond: Bond) -> usize {
        smiles
            .bond_matrix()
            .sparse_entries()
            .filter(|((row, column), entry)| *row < *column && entry.bond() == bond)
            .count()
    }

    fn count_aromatic_atoms(smiles: &Smiles) -> usize {
        smiles.nodes().iter().filter(|atom| atom.aromatic()).count()
    }

    fn aromatic_carbon() -> Atom {
        Atom::new_organic_subset(AtomSymbol::Element(Element::C), true)
    }

    fn smiles_from_edges(atom_nodes: Vec<Atom>, bond_edges: &[BondEdge]) -> Smiles {
        let mut builder = BondMatrixBuilder::with_capacity(bond_edges.len());
        for edge in bond_edges {
            builder.push_edge(edge.node_a(), edge.node_b(), edge.bond(), edge.ring_num()).unwrap();
        }
        let number_of_nodes = atom_nodes.len();
        Smiles::from_bond_matrix_parts(atom_nodes, builder.finish(number_of_nodes))
    }

    #[test]
    fn ringless_graph_returns_unchanged() {
        let smiles = Smiles::from_str("CC").expect("valid smiles");
        let kekulized = smiles.kekulize().expect("ringless graph should succeed");
        assert_eq!(kekulized, smiles);
    }

    #[test]
    fn benzene_kekulizes_to_explicit_bond_orders() {
        let smiles = Smiles::from_str("c1ccccc1").expect("valid aromatic benzene");
        let kekulized = smiles.kekulize().expect("benzene should kekulize");

        assert_eq!(count_bonds(&kekulized, Bond::Aromatic), 0);
        assert_eq!(count_bonds(&kekulized, Bond::Double), 3);
        assert_eq!(count_bonds(&kekulized, Bond::Single), 3);
        assert_eq!(count_aromatic_atoms(&kekulized), 0);
    }

    #[test]
    fn biphenyl_kekulizes_each_aromatic_component() {
        let smiles = Smiles::from_str("c1ccccc1-c2ccccc2").expect("valid aromatic biphenyl");
        let kekulized = smiles.kekulize().expect("biphenyl should kekulize");

        assert_eq!(count_bonds(&kekulized, Bond::Aromatic), 0);
        assert_eq!(count_bonds(&kekulized, Bond::Double), 6);
        assert_eq!(count_aromatic_atoms(&kekulized), 0);
    }

    #[test]
    fn aromatic_pyrrole_kekulizes_with_two_double_bonds() {
        let smiles = Smiles::from_str("[nH]1cccc1").expect("valid aromatic pyrrole");
        let kekulized = smiles.kekulize().expect("pyrrole should kekulize");

        assert_eq!(count_bonds(&kekulized, Bond::Aromatic), 0);
        assert_eq!(count_bonds(&kekulized, Bond::Double), 2);
        assert_eq!(count_bonds(&kekulized, Bond::Single), 3);
        assert_eq!(count_aromatic_atoms(&kekulized), 0);
    }

    #[test]
    fn pyridine_kekulizes_with_three_double_bonds() {
        let smiles = Smiles::from_str("n1ccccc1").expect("valid aromatic pyridine");
        let kekulized = smiles.kekulize().expect("pyridine should kekulize");

        assert_eq!(count_bonds(&kekulized, Bond::Aromatic), 0);
        assert_eq!(count_bonds(&kekulized, Bond::Double), 3);
        assert_eq!(count_bonds(&kekulized, Bond::Single), 3);
        assert_eq!(count_aromatic_atoms(&kekulized), 0);
    }

    #[test]
    fn aromatic_furan_kekulizes_with_two_double_bonds() {
        let smiles = Smiles::from_str("o1cccc1").expect("valid aromatic furan");
        let kekulized = smiles.kekulize().expect("furan should kekulize");

        assert_eq!(count_bonds(&kekulized, Bond::Aromatic), 0);
        assert_eq!(count_bonds(&kekulized, Bond::Double), 2);
        assert_eq!(count_bonds(&kekulized, Bond::Single), 3);
        assert_eq!(count_aromatic_atoms(&kekulized), 0);
    }

    #[test]
    fn aromatic_imidazole_kekulizes() {
        let smiles = Smiles::from_str("c1ncc[nH]1").expect("valid aromatic imidazole");
        let kekulized = smiles.kekulize().expect("imidazole should kekulize");

        assert_eq!(count_bonds(&kekulized, Bond::Aromatic), 0);
        assert_eq!(count_bonds(&kekulized, Bond::Double), 2);
        assert_eq!(count_bonds(&kekulized, Bond::Single), 3);
        assert_eq!(count_aromatic_atoms(&kekulized), 0);
    }

    #[test]
    fn aromatic_phenyl_anion_kekulizes() {
        let smiles = Smiles::from_str("c1cc[c-]cc1").expect("valid aromatic phenyl anion");
        let kekulized = smiles.kekulize().expect("phenyl anion should kekulize");

        assert_eq!(count_bonds(&kekulized, Bond::Aromatic), 0);
        assert_eq!(count_bonds(&kekulized, Bond::Double), 3);
        assert_eq!(count_aromatic_atoms(&kekulized), 0);
    }

    #[test]
    fn perceived_imidazole_retains_implicit_hydrogens_for_kekulization() {
        let smiles = Smiles::from_str("C1=CN=CN1").expect("valid Kekule imidazole");
        let expected_total_hydrogens = smiles
            .nodes()
            .iter()
            .enumerate()
            .map(|(atom_id, atom)| {
                atom.hydrogen_count() + smiles.implicit_hydrogen_count(atom_id).unwrap_or(0)
            })
            .collect::<Vec<_>>();
        let aromaticized = smiles
            .perceive_aromaticity()
            .expect("imidazole should aromaticize")
            .into_aromaticized();

        let aromatic_total_hydrogens = aromaticized
            .nodes()
            .iter()
            .enumerate()
            .map(|(atom_id, atom)| {
                atom.hydrogen_count() + aromaticized.implicit_hydrogen_count(atom_id).unwrap_or(0)
            })
            .collect::<Vec<_>>();

        assert_eq!(aromatic_total_hydrogens, expected_total_hydrogens);
        assert!(aromaticized.kekulize().is_ok());
    }

    #[test]
    fn preserve_source_mode_restores_original_kekule_graph() {
        let smiles = Smiles::from_str("C1=CN=CN1").expect("valid Kekule imidazole");
        let aromaticized = smiles
            .perceive_aromaticity()
            .expect("imidazole should aromaticize")
            .into_aromaticized();

        let restored =
            aromaticized.kekulize().expect("preserve-source kekulization should succeed");
        assert_eq!(restored, smiles);
    }

    #[test]
    fn standalone_kekulization_matches_explicit_mode() {
        let smiles = Smiles::from_str("c1ccccc1").expect("valid aromatic benzene");

        assert_eq!(
            smiles.kekulize_standalone().expect("standalone kekulization should succeed"),
            smiles
                .kekulize_with(KekulizationMode::Standalone)
                .expect("explicit standalone kekulization should succeed")
        );
    }

    #[test]
    fn perceived_indole_retains_implicit_hydrogens_for_kekulization() {
        let smiles = Smiles::from_str("C1=CC=C2C(=C1)C=CN2").expect("valid Kekule indole");
        let expected_total_hydrogens = smiles
            .nodes()
            .iter()
            .enumerate()
            .map(|(atom_id, atom)| {
                atom.hydrogen_count() + smiles.implicit_hydrogen_count(atom_id).unwrap_or(0)
            })
            .collect::<Vec<_>>();
        let aromaticized =
            smiles.perceive_aromaticity().expect("indole should aromaticize").into_aromaticized();

        let aromatic_total_hydrogens = aromaticized
            .nodes()
            .iter()
            .enumerate()
            .map(|(atom_id, atom)| {
                atom.hydrogen_count() + aromaticized.implicit_hydrogen_count(atom_id).unwrap_or(0)
            })
            .collect::<Vec<_>>();

        assert_eq!(aromatic_total_hydrogens, expected_total_hydrogens);
        assert!(aromaticized.kekulize().is_ok());
    }

    #[test]
    fn candidate_atoms_exclude_pyrrolic_heteroatoms() {
        let pyrrole = Smiles::from_str("[nH]1cccc1").expect("valid aromatic pyrrole");
        assert_eq!(candidate_atom_ids(&pyrrole), vec![1, 2, 3, 4]);

        let imidazole = Smiles::from_str("c1ncc[nH]1").expect("valid aromatic imidazole");
        assert_eq!(candidate_atom_ids(&imidazole), vec![0, 1, 2, 3]);
    }

    #[test]
    fn azulene_kekulizes() {
        let smiles = Smiles::from_str("c1cccc2cccc2c1").expect("valid aromatic azulene");
        let kekulized = smiles.kekulize().expect("azulene should kekulize");

        assert_eq!(count_bonds(&kekulized, Bond::Aromatic), 0);
        assert_eq!(count_bonds(&kekulized, Bond::Double), 5);
        assert_eq!(count_bonds(&kekulized, Bond::Single), 6);
        assert_eq!(count_aromatic_atoms(&kekulized), 0);
    }

    #[test]
    fn pyrrolic_nitrogen_without_hydrogen_fails_kekulization() {
        let smiles = Smiles::from_str("n1cccc1").expect("valid aromatic five-membered ring");

        assert_eq!(
            smiles.kekulize(),
            Err(KekulizationError::NoPerfectMatching { candidate_atom_count: 5 })
        );
    }

    #[test]
    fn odd_cycle_candidate_graph_has_no_perfect_matching() {
        let smiles = smiles_from_edges(
            vec![aromatic_carbon(), aromatic_carbon(), aromatic_carbon()],
            &[
                BondEdge::new(0, 1, Bond::Aromatic, None),
                BondEdge::new(1, 2, Bond::Aromatic, None),
                BondEdge::new(2, 0, Bond::Aromatic, None),
            ],
        );
        let aromatic_bonds = smiles
            .bond_matrix()
            .sparse_entries()
            .filter_map(|((row, column), entry)| {
                (row < column && entry.bond() == Bond::Aromatic).then_some(AromaticBond {
                    node_a: row,
                    node_b: column,
                    order: entry.order(),
                })
            })
            .collect::<Vec<_>>();
        let candidate_atom_ids = candidate_atom_ids(&smiles);
        let candidate_graph =
            KekulizationCandidateGraph::new(&smiles, &aromatic_bonds, &candidate_atom_ids);

        assert_eq!(
            candidate_graph.matched_original_orders(),
            Err(KekulizationError::NoPerfectMatching { candidate_atom_count: 3 })
        );
    }

    fn assert_roundtrip_for_policy(smiles: &str, policy: crate::smiles::AromaticityPolicy) {
        let smiles = Smiles::from_str(smiles).expect("valid smiles");
        let perceived =
            smiles.perceive_aromaticity_for(policy).expect("initial perception should succeed");
        let expected_atom_ids = perceived.assignment().atom_ids().to_vec();
        let expected_bond_edges = perceived.assignment().bond_edges().to_vec();
        let candidate_atom_ids = candidate_atom_ids(perceived.aromaticized());
        let atom_state = perceived
            .aromaticized()
            .nodes()
            .iter()
            .enumerate()
            .map(|(atom_id, atom)| {
                format!(
                    "{atom_id}: element={:?} aromatic={} charge={} explicit_h={} implicit_h={} degree={}",
                    atom.element(),
                    atom.aromatic(),
                    atom.charge_value(),
                    atom.hydrogen_count(),
                    perceived.aromaticized().implicit_hydrogen_count(atom_id).unwrap_or(0),
                    perceived.aromaticized().edges_for_node(atom_id).len(),
                )
            })
            .collect::<Vec<_>>();

        let kekulized = perceived.aromaticized().kekulize().unwrap_or_else(|error| {
            panic!(
                "kekulization should succeed: {error}; candidates={candidate_atom_ids:?}; atom_state={atom_state:?}"
            )
        });
        let reperceived =
            kekulized.perceive_aromaticity_for(policy).expect("reperception should succeed");

        assert_eq!(reperceived.assignment().atom_ids(), expected_atom_ids);
        assert_eq!(reperceived.assignment().bond_edges(), expected_bond_edges);
    }

    fn assert_rendered_aromatic_roundtrip_for_policy(
        smiles: &str,
        policy: crate::smiles::AromaticityPolicy,
    ) {
        let smiles = Smiles::from_str(smiles).expect("valid smiles");
        let perceived =
            smiles.perceive_aromaticity_for(policy).expect("initial perception should succeed");
        let expected_atom_ids = perceived.assignment().atom_ids().to_vec();
        let expected_bond_edges = perceived.assignment().bond_edges().to_vec();
        let aromatic_smiles = perceived.aromaticized().to_string();
        let reparsed = Smiles::from_str(&aromatic_smiles)
            .unwrap_or_else(|error| panic!("rendered aromatic SMILES should parse: {error}"));
        let candidate_atom_ids = candidate_atom_ids(&reparsed);
        let atom_state = reparsed
            .nodes()
            .iter()
            .enumerate()
            .map(|(atom_id, atom)| {
                format!(
                    "{atom_id}: element={:?} aromatic={} charge={} explicit_h={} implicit_h={} degree={}",
                    atom.element(),
                    atom.aromatic(),
                    atom.charge_value(),
                    atom.hydrogen_count(),
                    reparsed.implicit_hydrogen_count(atom_id).unwrap_or(0),
                    reparsed.edges_for_node(atom_id).len(),
                )
            })
            .collect::<Vec<_>>();

        let kekulized = reparsed.kekulize().unwrap_or_else(|error| {
            panic!(
                "standalone aromatic kekulization should succeed: {error}; aromatic_smiles={aromatic_smiles}; candidates={candidate_atom_ids:?}; atom_state={atom_state:?}"
            )
        });
        let reperceived =
            kekulized.perceive_aromaticity_for(policy).expect("reperception should succeed");

        assert_eq!(reperceived.assignment().atom_ids(), expected_atom_ids);
        assert_eq!(reperceived.assignment().bond_edges(), expected_bond_edges);
    }

    #[test]
    fn phosphorus_cage_roundtrips_under_default_perception() {
        assert_roundtrip_for_policy(
            "CP1P2P(P3P1P3P2C)C",
            crate::smiles::AromaticityPolicy::RdkitDefault,
        );
    }

    #[test]
    fn cyclopentadienyl_magnesium_complex_roundtrips_under_mdl_perception() {
        assert_roundtrip_for_policy(
            "C1=CC=[C-]C=C1.C1=CC=[C-]C=C1.[Mg+2]",
            crate::smiles::AromaticityPolicy::RdkitMdl,
        );
    }

    #[test]
    fn cyclopentadienyl_magnesium_complex_roundtrips_under_simple_perception() {
        assert_roundtrip_for_policy(
            "C1=CC=[C-]C=C1.C1=CC=[C-]C=C1.[Mg+2]",
            crate::smiles::AromaticityPolicy::RdkitSimple,
        );
    }

    #[test]
    fn rendered_aromatic_imidazole_roundtrips_under_default_perception() {
        assert_rendered_aromatic_roundtrip_for_policy(
            "C1=CN=CN1",
            crate::smiles::AromaticityPolicy::RdkitDefault,
        );
    }

    #[test]
    fn rendered_aromatic_phosphorus_cage_requires_preserved_source_under_default_perception() {
        let smiles = Smiles::from_str("CP1P2P(P3P1P3P2C)C").expect("valid smiles");
        let perceived = smiles
            .perceive_aromaticity_for(crate::smiles::AromaticityPolicy::RdkitDefault)
            .expect("initial perception should succeed");
        let aromatic_smiles = perceived.aromaticized().to_string();
        let reparsed =
            Smiles::from_str(&aromatic_smiles).expect("rendered aromatic smiles should parse");

        assert_eq!(
            reparsed.kekulize_standalone(),
            Err(KekulizationError::NoPerfectMatching { candidate_atom_count: 0 })
        );
    }

    #[test]
    fn rendered_aromatic_cyclopentadienyl_magnesium_complex_roundtrips_under_mdl_perception() {
        assert_rendered_aromatic_roundtrip_for_policy(
            "C1=CC=[C-]C=C1.C1=CC=[C-]C=C1.[Mg+2]",
            crate::smiles::AromaticityPolicy::RdkitMdl,
        );
    }

    #[test]
    fn rendered_aromatic_cyclopentadienyl_magnesium_complex_roundtrips_under_simple_perception() {
        assert_rendered_aromatic_roundtrip_for_policy(
            "C1=CC=[C-]C=C1.C1=CC=[C-]C=C1.[Mg+2]",
            crate::smiles::AromaticityPolicy::RdkitSimple,
        );
    }

    #[test]
    fn fuzz_regression_star_hash_2o_star_star_2() {
        let smiles = Smiles::from_str("*#2O**2").expect("fuzz input should parse");

        for policy in [
            crate::smiles::AromaticityPolicy::RdkitDefault,
            crate::smiles::AromaticityPolicy::RdkitSimple,
            crate::smiles::AromaticityPolicy::RdkitMdl,
        ] {
            let Ok(perception) = smiles.perceive_aromaticity_for(policy) else {
                continue;
            };

            let preserved = perception
                .kekulize()
                .expect("preserve-source kekulization should succeed after perception");
            assert_eq!(&preserved, &smiles);

            if perception.status() == crate::smiles::AromaticityStatus::Unsupported {
                continue;
            }

            if perception.assignment().bond_edges().is_empty() {
                assert_eq!(
                    perception
                        .kekulize_standalone()
                        .expect("non-aromatic standalone kekulization should succeed"),
                    *perception.aromaticized(),
                    "policy={policy:?}; input={smiles}; aromaticized={}",
                    perception.aromaticized()
                );
            } else {
                match perception.kekulize_standalone() {
                    Ok(standalone) => {
                        let reperceived = standalone
                            .perceive_aromaticity_for(policy)
                            .expect("reperception should succeed after standalone kekulization");

                        assert_eq!(
                            reperceived.assignment().atom_ids(),
                            perception.assignment().atom_ids(),
                            "policy={policy:?}; input={smiles}; aromaticized={}; standalone={}",
                            perception.aromaticized(),
                            standalone
                        );
                        assert_eq!(
                            reperceived.assignment().bond_edges(),
                            perception.assignment().bond_edges(),
                            "policy={policy:?}; input={smiles}; aromaticized={}; standalone={}",
                            perception.aromaticized(),
                            standalone
                        );
                    }
                    Err(error) => {
                        assert_eq!(
                            error,
                            KekulizationError::NoPerfectMatching { candidate_atom_count: 0 },
                            "policy={policy:?}; input={smiles}; aromaticized={}",
                            perception.aromaticized()
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn fuzz_regression_complex_wildcard_ring_case() {
        let smiles = Smiles::from_str(
            "*13**#C*#O*#O*12*Br**13**-N5NNNNN5NNN5NNNNN5NNN*N5N=NNNN5N*N5N=NNNN5N1*2#2O**2*",
        )
        .expect("fuzz input should parse");

        for policy in [
            crate::smiles::AromaticityPolicy::RdkitDefault,
            crate::smiles::AromaticityPolicy::RdkitSimple,
            crate::smiles::AromaticityPolicy::RdkitMdl,
        ] {
            let Ok(perception) = smiles.perceive_aromaticity_for(policy) else {
                continue;
            };

            let preserved = perception
                .kekulize()
                .expect("preserve-source kekulization should succeed after perception");
            assert_eq!(&preserved, &smiles);

            if perception.status() == crate::smiles::AromaticityStatus::Unsupported {
                continue;
            }

            if perception.assignment().bond_edges().is_empty() {
                assert_eq!(
                    perception
                        .kekulize_standalone()
                        .expect("non-aromatic standalone kekulization should succeed"),
                    *perception.aromaticized(),
                    "policy={policy:?}; input={smiles}; aromaticized={}",
                    perception.aromaticized()
                );
            } else {
                match perception.kekulize_standalone() {
                    Ok(standalone) => {
                        let reperceived = standalone
                            .perceive_aromaticity_for(policy)
                            .expect("reperception should succeed after standalone kekulization");

                        assert_eq!(
                            reperceived.assignment().atom_ids(),
                            perception.assignment().atom_ids(),
                            "policy={policy:?}; input={smiles}; aromaticized={}; standalone={}",
                            perception.aromaticized(),
                            standalone
                        );
                        assert_eq!(
                            reperceived.assignment().bond_edges(),
                            perception.assignment().bond_edges(),
                            "policy={policy:?}; input={smiles}; aromaticized={}; standalone={}",
                            perception.aromaticized(),
                            standalone
                        );
                    }
                    Err(error) => {
                        assert_eq!(
                            error,
                            KekulizationError::NoPerfectMatching { candidate_atom_count: 0 },
                            "policy={policy:?}; input={smiles}; aromaticized={}",
                            perception.aromaticized()
                        );
                    }
                }
            }
        }
    }
}
