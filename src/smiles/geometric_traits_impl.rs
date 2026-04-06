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
fn is_row_major_sorted(entries: &[PendingBond]) -> bool {
    entries.windows(2).all(|window| window[0].row_major_key() <= window[1].row_major_key())
}

impl Smiles {
    #[inline]
    #[must_use]
    pub(crate) fn from_bond_matrix_parts(atom_nodes: Vec<Atom>, bond_matrix: BondMatrix) -> Self {
        Self { atom_nodes, bond_matrix }
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

    use geometric_traits::traits::{GraphSimilarities, McesBuilder};

    use super::*;
    use crate::bond::ring_num::RingNum;

    fn assert_similarity_close(actual: impl Into<f64>, expected: f64) {
        let actual = actual.into();
        let difference = (actual - expected).abs();

        assert!(
            difference < 1.0e-12,
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
}
