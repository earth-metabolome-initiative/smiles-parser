use alloc::boxed::Box;
use alloc::vec::Vec;
use core::hash::{Hash, Hasher};

use geometric_traits::{
    impls::{SymmetricCSR2D, ValuedCSR2D},
    traits::{Graph, MonopartiteGraph, MonoplexGraph, SizedSparseMatrix},
};
use hashbrown::HashSet;

use super::Smiles;
use crate::{
    atom::Atom,
    bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
    errors::SmilesError,
};

/// The symmetric valued sparse matrix storing SMILES bonds.
pub type BondMatrix = SymmetricCSR2D<ValuedCSR2D<usize, usize, usize, BondEntry>>;

#[derive(Clone, Copy, Debug)]
/// Value stored for each bond in the symmetric adjacency matrix.
pub struct BondEntry {
    bond: Bond,
    ring_num: Option<RingNum>,
    order: usize,
}

impl BondEntry {
    /// Creates a new stored bond value for the adjacency matrix.
    #[inline]
    #[must_use]
    pub const fn new(bond: Bond, ring_num: Option<RingNum>, order: usize) -> Self {
        Self { bond, ring_num, order }
    }

    /// Returns the bond type stored for this adjacency entry.
    #[inline]
    #[must_use]
    pub const fn bond(self) -> Bond {
        self.bond
    }

    /// Returns the ring number stored for this adjacency entry, if any.
    #[inline]
    #[must_use]
    pub const fn ring_num(self) -> Option<RingNum> {
        self.ring_num
    }

    #[inline]
    #[must_use]
    pub(crate) const fn order(self) -> usize {
        self.order
    }

    #[inline]
    #[must_use]
    pub(crate) fn to_bond_edge(self, node_a: usize, node_b: usize) -> BondEdge {
        BondEdge::new(node_a, node_b, self.bond, self.ring_num)
    }

    #[inline]
    #[must_use]
    pub(crate) const fn with_bond(mut self, bond: Bond) -> Self {
        self.bond = bond;
        self
    }

    #[inline]
    #[must_use]
    pub(crate) const fn with_order(mut self, order: usize) -> Self {
        self.order = order;
        self
    }
}

impl PartialEq for BondEntry {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.bond == other.bond
    }
}

impl Eq for BondEntry {}

impl Hash for BondEntry {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bond.hash(state);
    }
}

#[derive(Debug, Default)]
pub(crate) struct BondMatrixBuilder {
    entries: Vec<PendingBond>,
    seen_edges: HashSet<(usize, usize)>,
}

impl BondMatrixBuilder {
    #[inline]
    #[must_use]
    pub(crate) fn with_capacity(number_of_edges: usize) -> Self {
        Self {
            entries: Vec::with_capacity(number_of_edges),
            seen_edges: HashSet::with_capacity(number_of_edges),
        }
    }

    #[inline]
    #[must_use]
    pub(crate) fn contains_edge(&self, node_a: usize, node_b: usize) -> bool {
        let (row, column) = Smiles::edge_key(node_a, node_b);
        self.seen_edges.contains(&(row, column))
    }

    #[inline]
    pub(crate) fn push_edge(
        &mut self,
        node_a: usize,
        node_b: usize,
        bond: Bond,
        ring_num: Option<RingNum>,
    ) -> Result<(), SmilesError> {
        let (row, column) = Smiles::edge_key(node_a, node_b);
        if row == column {
            return Err(SmilesError::SelfLoopEdge(row));
        }
        if !self.seen_edges.insert((row, column)) {
            return Err(SmilesError::DuplicateEdge(row, column));
        }

        let order = self.entries.len();
        self.entries.push(PendingBond::new(row, column, BondEntry::new(bond, ring_num, order)));
        Ok(())
    }

    #[inline]
    #[must_use]
    pub(crate) fn finish(self, number_of_nodes: usize) -> BondMatrix {
        build_bond_matrix(number_of_nodes, self.entries)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct PendingBond {
    row: usize,
    column: usize,
    entry: BondEntry,
}

impl PendingBond {
    #[inline]
    #[must_use]
    const fn new(row: usize, column: usize, entry: BondEntry) -> Self {
        Self { row, column, entry }
    }

    #[inline]
    #[must_use]
    const fn row_major_key(self) -> (usize, usize) {
        (self.row, self.column)
    }

    #[inline]
    #[must_use]
    const fn into_entry(self) -> (usize, usize, BondEntry) {
        (self.row, self.column, self.entry)
    }
}

#[inline]
#[must_use]
fn build_bond_matrix(number_of_nodes: usize, mut entries: Vec<PendingBond>) -> BondMatrix {
    reassign_rdkit_bond_orders(&mut entries);
    if !is_row_major_sorted(&entries) {
        entries.sort_unstable_by_key(|bond| bond.row_major_key());
    }
    BondMatrix::from_sorted_upper_triangular_entries(
        number_of_nodes,
        entries.into_iter().map(PendingBond::into_entry),
    )
    .unwrap_or_else(|_| {
        unreachable!("bond entries are unique, upper-triangular, and row-major sorted")
    })
}

#[inline]
fn reassign_rdkit_bond_orders(entries: &mut [PendingBond]) {
    let mut reordered = entries.iter_mut().collect::<Vec<_>>();
    reordered.sort_unstable_by_key(|pending| {
        (pending.entry.ring_num().is_some(), pending.entry.order())
    });
    for (order, pending) in reordered.into_iter().enumerate() {
        pending.entry = pending.entry.with_order(order);
    }
}

#[inline]
fn is_row_major_sorted(entries: &[PendingBond]) -> bool {
    entries.windows(2).all(|window| window[0].row_major_key() <= window[1].row_major_key())
}

impl Smiles {
    #[inline]
    #[must_use]
    pub(crate) fn from_bond_matrix_parts(atom_nodes: Vec<Atom>, bond_matrix: BondMatrix) -> Self {
        Self::from_bond_matrix_parts_with_caches(atom_nodes, bond_matrix, None, None)
    }

    #[inline]
    #[must_use]
    pub(crate) fn from_bond_matrix_parts_with_implicit_hydrogen_cache(
        atom_nodes: Vec<Atom>,
        bond_matrix: BondMatrix,
        implicit_hydrogen_cache: Vec<u8>,
    ) -> Self {
        Self::from_bond_matrix_parts_with_caches(
            atom_nodes,
            bond_matrix,
            Some(implicit_hydrogen_cache),
            None,
        )
    }

    #[inline]
    #[must_use]
    pub(crate) fn from_bond_matrix_parts_with_caches(
        atom_nodes: Vec<Atom>,
        bond_matrix: BondMatrix,
        implicit_hydrogen_cache: Option<Vec<u8>>,
        kekulization_source: Option<Box<Self>>,
    ) -> Self {
        Self { atom_nodes, bond_matrix, implicit_hydrogen_cache, kekulization_source }
    }

    /// Returns the symmetric valued sparse matrix storing the graph bonds.
    #[inline]
    #[must_use]
    pub fn bond_matrix(&self) -> &BondMatrix {
        &self.bond_matrix
    }

    /// Returns the number of unique chemical bonds in the graph.
    #[inline]
    #[must_use]
    pub fn number_of_bonds(&self) -> usize {
        self.bond_matrix.number_of_defined_values() / 2
    }
}

impl Graph for Smiles {
    #[inline]
    fn has_nodes(&self) -> bool {
        !self.atom_nodes.is_empty()
    }

    #[inline]
    fn has_edges(&self) -> bool {
        self.number_of_bonds() > 0
    }
}

impl MonoplexGraph for Smiles {
    type Edge = <BondMatrix as geometric_traits::traits::Edges>::Edge;
    type Edges = BondMatrix;

    #[inline]
    fn edges(&self) -> &Self::Edges {
        &self.bond_matrix
    }
}

impl MonopartiteGraph for Smiles {
    type NodeId = usize;
    type NodeSymbol = Atom;
    type Nodes = Vec<Atom>;

    #[inline]
    fn nodes_vocabulary(&self) -> &Self::Nodes {
        &self.atom_nodes
    }
}

#[cfg(test)]
mod tests {
    use core::str::FromStr;
    use std::{
        collections::hash_map::DefaultHasher,
        hash::{Hash, Hasher},
    };

    use geometric_traits::traits::{Graph, GraphSimilarities, McesBuilder};

    use super::*;
    use crate::bond::ring_num::RingNum;

    fn assert_similarity_close(actual: impl Into<f64>, expected: f64) {
        let actual = actual.into();
        let difference = (actual - expected).abs();

        assert!(
            difference < 1.0e-6,
            "expected Johnson similarity {expected}, got {actual} (diff {difference})"
        );
    }

    #[test]
    fn bond_entry_equality_ignores_ring_digits_and_uses_bond_kind() {
        let first = BondEntry::new(Bond::Double, Some(RingNum::try_new(1).unwrap()), 0);
        let second = BondEntry::new(Bond::Double, Some(RingNum::try_new(9).unwrap()), 17);
        let third = BondEntry::new(Bond::Double, None, 99);
        let fourth = BondEntry::new(Bond::Single, Some(RingNum::try_new(1).unwrap()), 0);

        assert_eq!(first, second);
        assert_eq!(first, third);
        assert_ne!(first, fourth);
    }

    #[test]
    fn bond_entry_hash_ignores_ring_digits_like_equality() {
        let first = BondEntry::new(Bond::Double, Some(RingNum::try_new(1).unwrap()), 0);
        let second = BondEntry::new(Bond::Double, Some(RingNum::try_new(9).unwrap()), 17);
        let mut first_hasher = DefaultHasher::new();
        let mut second_hasher = DefaultHasher::new();

        first.hash(&mut first_hasher);
        second.hash(&mut second_hasher);

        assert_eq!(first_hasher.finish(), second_hasher.finish());
    }

    #[test]
    fn bond_matrix_builder_rejects_duplicate_edges() {
        let mut builder = BondMatrixBuilder::with_capacity(2);

        builder.push_edge(0, 1, Bond::Single, None).unwrap();
        let error = builder.push_edge(1, 0, Bond::Double, None).unwrap_err();

        assert_eq!(error, SmilesError::DuplicateEdge(0, 1));
    }

    #[test]
    fn empty_smiles_reports_no_nodes_or_edges() {
        let smiles = Smiles::default();

        assert_eq!(smiles.number_of_bonds(), 0);
        assert!(!Graph::has_nodes(&smiles));
        assert!(!Graph::has_edges(&smiles));
    }

    struct DirectionalParityCase {
        name: &'static str,
        smiles1: &'static str,
        smiles2: &'static str,
        raw_edges: usize,
        raw_similarity: f64,
        collapsed_edges: usize,
        collapsed_similarity: f64,
    }

    #[test]
    fn labeled_mces_distinguishes_atom_identity() {
        let ethanol = Smiles::from_str("CCO").unwrap();
        let ethylamine = Smiles::from_str("CCN").unwrap();
        let dimethyl_ether = Smiles::from_str("COC").unwrap();

        let result_one = McesBuilder::new(&ethanol, &ethylamine).compute_labeled();
        let result_two = McesBuilder::new(&ethanol, &dimethyl_ether).compute_labeled();

        assert_eq!(result_one.matched_edges().len(), 1);
        assert_eq!(result_two.matched_edges().len(), 1);
        assert!(result_one.johnson_similarity() < 1.0);
        assert!(result_two.johnson_similarity() < 1.0);
    }

    #[test]
    fn labeled_mces_ignores_bracket_and_explicit_hydrogen_spelling() {
        let implicit = Smiles::from_str("CC").unwrap();
        let explicit = Smiles::from_str("[CH3][CH3]").unwrap();

        let result = McesBuilder::new(&implicit, &explicit).compute_labeled();

        assert_eq!(result.matched_edges().len(), 1);
        assert_similarity_close(result.johnson_similarity(), 1.0);
    }

    #[test]
    fn labeled_mces_rejects_topology_only_ring_false_positive() {
        let benzene = Smiles::from_str("c1ccccc1").unwrap();
        let pyridine = Smiles::from_str("c1ccncc1").unwrap();

        let result = McesBuilder::new(&benzene, &pyridine).compute_labeled();

        assert_eq!(result.matched_edges().len(), 4);
        assert!(result.johnson_similarity() < 1.0);
    }

    #[test]
    fn labeled_mces_ignores_literal_ring_digits_in_bond_labels() {
        let ring_one = Smiles::from_str("C1CCCCC1").unwrap();
        let ring_two = Smiles::from_str("C2CCCCC2").unwrap();

        let result = McesBuilder::new(&ring_one, &ring_two).compute_labeled();

        assert_eq!(result.matched_edges().len(), 6);
        assert_similarity_close(result.johnson_similarity(), 1.0);
    }

    #[test]
    fn directional_bond_collapse_recovers_eight_massspecgym_directional_cases() {
        let cases = [
            DirectionalParityCase {
                name: "massspecgym_default_0012",
                smiles1: "CCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
                smiles2: "CCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC",
                raw_edges: 42,
                raw_similarity: 0.728_139,
                collapsed_edges: 46,
                collapsed_similarity: 0.868_206,
            },
            DirectionalParityCase {
                name: "massspecgym_default_0015",
                smiles1: "CCCCCCC=CCCCCCCCC(=O)O",
                smiles2: "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC",
                raw_edges: 15,
                raw_similarity: 0.723_588,
                collapsed_edges: 17,
                collapsed_similarity: 0.813_953,
            },
            DirectionalParityCase {
                name: "massspecgym_default_0051",
                smiles1: "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCC(=O)O",
                smiles2: "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
                raw_edges: 41,
                raw_similarity: 0.723_406,
                collapsed_edges: 42,
                collapsed_similarity: 0.758_689,
            },
            DirectionalParityCase {
                name: "massspecgym_default_0059",
                smiles1: "CCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC",
                smiles2: "CCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC/C=C\\CCCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C",
                raw_edges: 45,
                raw_similarity: 0.882_461,
                collapsed_edges: 47,
                collapsed_similarity: 0.960_004,
            },
            DirectionalParityCase {
                name: "massspecgym_default_0069",
                smiles1: "CCCCCC/C=C\\CCCCCCCC(=O)O",
                smiles2: "CCCCCCC=CCCCCCCCCC(=O)NCC(=O)O",
                raw_edges: 14,
                raw_similarity: 0.650_159,
                collapsed_edges: 16,
                collapsed_similarity: 0.733_968,
            },
            DirectionalParityCase {
                name: "massspecgym_default_0090",
                smiles1: "CCCCCCCCCCCC(=O)OC(CCCCCC)C/C=C/CCCCCCCC(=O)O",
                smiles2: "CCCCCCCCCCCCCCCC(=O)OC(CCCCCCC)CCCCCCCCCCC(=O)O",
                raw_edges: 30,
                raw_similarity: 0.745_106,
                collapsed_edges: 32,
                collapsed_similarity: 0.844_350,
            },
            DirectionalParityCase {
                name: "massspecgym_default_0097",
                smiles1: "CCCCCCCC/C=C\\CCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC/C=C\\CCCCCCCC)COP(=O)(O)OCCN",
                smiles2: "CCCCCCCCCCCCCCCC/C=C\\OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC",
                raw_edges: 42,
                raw_similarity: 0.720_961,
                collapsed_edges: 46,
                collapsed_similarity: 0.854_829,
            },
            DirectionalParityCase {
                name: "massspecgym_default_0098",
                smiles1: "CCCCC/C=C\\CCCCCCCCCCCCC(=O)O",
                smiles2: "CCCCCC=CCC=CCCCCCCCCC(=O)OCCCCCCCCCCCCCCC(=O)O",
                raw_edges: 19,
                raw_similarity: 0.521_240,
                collapsed_edges: 20,
                collapsed_similarity: 0.546_977,
            },
        ];

        for case in cases {
            let smiles1 = Smiles::from_str(case.smiles1).unwrap();
            let smiles2 = Smiles::from_str(case.smiles2).unwrap();
            let raw = McesBuilder::new(&smiles1, &smiles2).compute_labeled();
            assert_eq!(raw.matched_edges().len(), case.raw_edges, "{}", case.name);
            assert_similarity_close(raw.johnson_similarity(), case.raw_similarity);

            let collapsed1 = smiles1.with_directional_bonds_collapsed();
            let collapsed2 = smiles2.with_directional_bonds_collapsed();
            let collapsed = McesBuilder::new(&collapsed1, &collapsed2).compute_labeled();
            assert_eq!(collapsed.matched_edges().len(), case.collapsed_edges, "{}", case.name);
            assert_similarity_close(collapsed.johnson_similarity(), case.collapsed_similarity);
        }
    }

    #[test]
    fn directional_bond_collapse_improves_but_does_not_fully_resolve_massspecgym_0023() {
        let smiles1 = Smiles::from_str(
            "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCC/C=C\\CCCCCC",
        )
        .unwrap();
        let smiles2 = Smiles::from_str(
            "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
        )
        .unwrap();

        let raw = McesBuilder::new(&smiles1, &smiles2).compute_labeled();
        let collapsed1 = smiles1.with_directional_bonds_collapsed();
        let collapsed2 = smiles2.with_directional_bonds_collapsed();
        let collapsed = McesBuilder::new(&collapsed1, &collapsed2).compute_labeled();

        assert_eq!(raw.matched_edges().len(), 43);
        assert_similarity_close(raw.johnson_similarity(), 0.780_422);
        assert_eq!(collapsed.matched_edges().len(), 44);
        assert_similarity_close(collapsed.johnson_similarity(), 0.797_861_065_613_257_5);
        assert!(collapsed.matched_edges().len() > raw.matched_edges().len());
        assert!(collapsed.johnson_similarity() > raw.johnson_similarity());
    }
}
