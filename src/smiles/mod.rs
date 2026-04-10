//! Represents a parsed SMILES graph.
//!
//! A [`Smiles`] value stores atoms as [`Atom`] values and bonds in a
//! symmetric valued sparse matrix.
//!
//! # Examples
//!
//! SMILES strings can be parsed into a graph, inspected, and rendered back
//! into a SMILES string.
//!
//! ```rust
//! use core::str::FromStr;
//!
//! use smiles_parser::prelude::Smiles;
//!
//! let source = "CC";
//! let smiles: Smiles = source.parse()?;
//!
//! assert_eq!(smiles.nodes().len(), 2);
//! assert_eq!(smiles.number_of_bonds(), 1);
//! assert_eq!(smiles.to_string(), "CC");
//!
//! # Ok::<(), smiles_parser::errors::SmilesErrorWithSpan>(())
//! ```
use alloc::{string::String, vec::Vec};
use core::fmt;

use geometric_traits::traits::{
    BiconnectedComponents, SizedSparseMatrix2D, SizedSparseValuedMatrixRef, SparseMatrix2D,
    SparseValuedMatrix2DRef, SparseValuedMatrixRef,
};

use crate::{
    atom::Atom,
    bond::bond_edge::BondEdge,
    errors::SmilesError,
    traversal::{render_visitor::RenderVisitor, walker::walk},
};

mod aromaticity;
mod from_str;
mod geometric_traits_impl;
mod implicit_hydrogens;
mod rdkit_symm_sssr;

pub(crate) use self::geometric_traits_impl::BondMatrixBuilder;
pub use self::{
    aromaticity::{
        AromaticityAssignment, AromaticityAssignmentApplicationError, AromaticityDiagnostic,
        AromaticityModel, AromaticityPerception, AromaticityPolicy, AromaticityRingFamilyKind,
        AromaticityStatus, RdkitDefaultAromaticity,
    },
    geometric_traits_impl::{BondEntry, BondMatrix},
};

/// Error raised while deriving ring membership from a [`Smiles`] graph.
pub type RingMembershipError = geometric_traits::errors::MonopartiteError<Smiles>;

/// Status flags describing whether symmetrized SSSR completed exactly or
/// required a fallback path.
#[derive(Debug, Copy, Clone, Default, PartialEq, Eq, Hash)]
pub struct SymmSssrStatus {
    used_fallback: bool,
    hit_queue_cutoff: bool,
}

impl SymmSssrStatus {
    /// Returns whether the SSSR search completed without fallback or cutoff.
    #[must_use]
    pub fn is_complete(self) -> bool {
        !self.used_fallback && !self.hit_queue_cutoff
    }

    /// Returns whether the search fell back to the approximate ring finder.
    #[must_use]
    pub fn used_fallback(self) -> bool {
        self.used_fallback
    }

    /// Returns whether the search hit the BFS queue cutoff.
    #[must_use]
    pub fn hit_queue_cutoff(self) -> bool {
        self.hit_queue_cutoff
    }
}

/// Canonicalized symmetrized SSSR cycles together with the search status.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SymmSssrResult {
    cycles: Vec<Vec<usize>>,
    status: SymmSssrStatus,
}

impl SymmSssrResult {
    /// Returns the canonicalized cycle set.
    #[must_use]
    pub fn cycles(&self) -> &[Vec<usize>] {
        &self.cycles
    }

    /// Returns the search status for the cycle set.
    #[must_use]
    pub fn status(&self) -> SymmSssrStatus {
        self.status
    }
}

/// Canonicalized ring-membership summary for a [`Smiles`] graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingMembership {
    atom_ids: Vec<usize>,
    bond_edges: Vec<[usize; 2]>,
}

impl RingMembership {
    /// Returns the atom ids that belong to at least one ring.
    #[inline]
    #[must_use]
    pub fn atom_ids(&self) -> &[usize] {
        &self.atom_ids
    }

    /// Returns the bond edges that belong to at least one ring.
    #[inline]
    #[must_use]
    pub fn bond_edges(&self) -> &[[usize; 2]] {
        &self.bond_edges
    }

    /// Returns whether the given atom id belongs to at least one ring.
    #[inline]
    #[must_use]
    pub fn contains_atom(&self, atom_id: usize) -> bool {
        self.atom_ids.binary_search(&atom_id).is_ok()
    }

    /// Returns whether the given edge belongs to at least one ring.
    #[inline]
    #[must_use]
    pub fn contains_edge(&self, node_a: usize, node_b: usize) -> bool {
        let edge = if node_a < node_b { [node_a, node_b] } else { [node_b, node_a] };
        self.bond_edges.binary_search(&edge).is_ok()
    }
}

/// Represents a parsed SMILES graph.
#[derive(Debug, PartialEq)]
pub struct Smiles {
    atom_nodes: Vec<Atom>,
    bond_matrix: BondMatrix,
}

impl Smiles {
    /// Creates a new empty [`Smiles`] graph.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self { atom_nodes: Vec::new(), bond_matrix: BondMatrix::default() }
    }

    /// Returns a slice of all parsed [`Atom`] values.
    #[inline]
    #[must_use]
    pub fn nodes(&self) -> &[Atom] {
        &self.atom_nodes
    }

    /// Returns the atom with the given positional index, if present.
    #[inline]
    #[must_use]
    pub fn node_by_id(&self, id: usize) -> Option<&Atom> {
        self.atom_nodes.get(id)
    }

    /// Returns a normalized edge key with node IDs in ascending order.
    #[inline]
    #[must_use]
    pub fn edge_key(node_a: usize, node_b: usize) -> (usize, usize) {
        if node_a < node_b { (node_a, node_b) } else { (node_b, node_a) }
    }

    /// Returns the bond connecting the given pair of node ids, if present.
    #[inline]
    #[must_use]
    pub fn edge_for_node_pair(&self, nodes: (usize, usize)) -> Option<BondEdge> {
        let (row, column) = Self::edge_key(nodes.0, nodes.1);
        let rank = self.bond_matrix.try_rank(row, column)?;
        let entry = *self.bond_matrix.select_value_ref(rank);
        Some(entry.to_bond_edge(row, column))
    }

    /// Returns the bonds incident to the provided node id.
    #[inline]
    #[must_use]
    pub fn edges_for_node(&self, id: usize) -> Vec<BondEdge> {
        if id >= self.atom_nodes.len() {
            return Vec::new();
        }

        self.bond_matrix
            .sparse_row(id)
            .zip(self.bond_matrix.sparse_row_values_ref(id))
            .map(|(other, entry)| entry.to_bond_edge(id, other))
            .collect()
    }

    /// Returns the atoms and bonds that belong to at least one ring.
    ///
    /// Ring membership is derived from the cyclic biconnected components of
    /// the graph, so chain attachments and bridges are excluded.
    ///
    /// # Errors
    /// Returns a [`RingMembershipError`] if the underlying graph violates the
    /// simple undirected assumptions of the biconnected-components algorithm.
    pub fn try_ring_membership(&self) -> Result<RingMembership, RingMembershipError> {
        let decomposition = self.biconnected_components()?;
        let mut atom_ids = Vec::new();
        let mut bond_edges = Vec::new();

        for component_id in decomposition.cyclic_biconnected_component_ids() {
            for edge in decomposition.edge_biconnected_component(component_id) {
                let bond_edge =
                    if edge[0] <= edge[1] { [edge[0], edge[1]] } else { [edge[1], edge[0]] };
                atom_ids.push(bond_edge[0]);
                atom_ids.push(bond_edge[1]);
                bond_edges.push(bond_edge);
            }
        }

        atom_ids.sort_unstable();
        atom_ids.dedup();
        bond_edges.sort_unstable();
        bond_edges.dedup();

        Ok(RingMembership { atom_ids, bond_edges })
    }

    /// Returns the atoms and bonds that belong to at least one ring.
    ///
    /// This is the infallible convenience wrapper around
    /// [`Smiles::try_ring_membership`]. The underlying biconnected-components
    /// algorithm only fails on unsupported graph shapes such as self-loops,
    /// which valid parsed [`Smiles`] graphs do not contain.
    #[inline]
    #[must_use]
    pub fn ring_membership(&self) -> RingMembership {
        self.try_ring_membership().unwrap_or_else(|_| {
            unreachable!("parsed SMILES graphs satisfy the ring-membership preconditions")
        })
    }

    /// Returns the canonicalized symmetric-SSSR-style cycle set together with
    /// search completeness status.
    ///
    /// The current implementation mirrors `RDKit`'s ring-finding pipeline,
    /// including its approximate fallback behavior on some highly fused
    /// systems.
    #[must_use]
    pub fn symm_sssr_result(&self) -> SymmSssrResult {
        rdkit_symm_sssr::symmetrize_sssr_with_status(self)
    }

    /// Returns a graph with directional single bonds collapsed to ordinary
    /// single bonds.
    ///
    /// This is intended for chemistry-facing matching paths where bond order
    /// should be compared without carrying slash/backslash directional notation
    /// through the MCES bond type.
    #[inline]
    #[must_use]
    pub fn with_directional_bonds_collapsed(&self) -> Self {
        let bond_matrix = BondMatrix::from_sorted_upper_triangular_entries(
            self.atom_nodes.len(),
            self.bond_matrix.sparse_entries().filter_map(|((row, column), entry)| {
                (row < column).then_some((
                    row,
                    column,
                    entry.with_bond(entry.bond().without_direction()),
                ))
            }),
        )
        .unwrap_or_else(|_| unreachable!("existing bond matrix entries are already valid"));

        Self { atom_nodes: self.atom_nodes.clone(), bond_matrix }
    }

    /// Renders the graph back into a valid SMILES string.
    ///
    /// # Errors
    /// - Returns a [`SmilesError`] if traversal fails.
    pub fn render(&self) -> Result<String, SmilesError> {
        self.render_visitor().map(RenderVisitor::into_string)
    }
    fn render_visitor(&self) -> Result<RenderVisitor, SmilesError> {
        let mut render_visitor =
            RenderVisitor::with_capacity(self.nodes().len(), self.number_of_bonds());
        walk(self, &mut render_visitor)?;
        Ok(render_visitor)
    }
}

fn canonicalize_cycle(cycle: &[usize]) -> Vec<usize> {
    let smallest_atom_id =
        *cycle.iter().min().unwrap_or_else(|| unreachable!("Johnson only yields non-empty cycles"));
    let start_index = cycle
        .iter()
        .position(|&atom_id| atom_id == smallest_atom_id)
        .unwrap_or_else(|| unreachable!("the smallest atom id is taken from the same cycle"));

    let forward =
        cycle[start_index..].iter().chain(cycle[..start_index].iter()).copied().collect::<Vec<_>>();

    let reversed = cycle.iter().rev().copied().collect::<Vec<_>>();
    let reversed_start_index = reversed
        .iter()
        .position(|&atom_id| atom_id == smallest_atom_id)
        .unwrap_or_else(|| unreachable!("the smallest atom id is taken from the same cycle"));
    let backward = reversed[reversed_start_index..]
        .iter()
        .chain(reversed[..reversed_start_index].iter())
        .copied()
        .collect::<Vec<_>>();

    forward.min(backward)
}

fn cycle_edges(cycle: &[usize]) -> Vec<[usize; 2]> {
    let number_of_nodes = cycle.len();
    (0..number_of_nodes)
        .map(|index| {
            let next_index = (index + 1) % number_of_nodes;
            let left = cycle[index];
            let right = cycle[next_index];
            if left <= right { [left, right] } else { [right, left] }
        })
        .collect()
}

#[cfg(test)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RingComponent {
    atom_ids: Vec<usize>,
    bond_edges: Vec<[usize; 2]>,
}

impl Default for Smiles {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Smiles {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let render_visitor = self.render_visitor().map_err(|_| fmt::Error)?;
        render_visitor.write_into_formatter(f)
    }
}

#[cfg(test)]
mod tests {
    use alloc::{string::ToString, vec::Vec};
    use std::str::FromStr;

    use elements_rs::Element;

    use super::{
        AromaticityAssignment, AromaticityAssignmentApplicationError, AromaticityDiagnostic,
        AromaticityModel, AromaticityPolicy, AromaticityStatus, BondMatrixBuilder, RingMembership,
        Smiles,
    };
    use crate::{
        atom::{Atom, atom_symbol::AtomSymbol},
        bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
        errors::SmilesError,
    };

    const LARGE_POLYCYCLE_FRONTIER_CASE: &str = "CC1=C2CC3=C(C4=C5C=C3[C@@H]6C2=CC7=C1CC8=C(C9=C1C=C8[C@@H]7CCC[C@@H]2C3=C7CC8=C2C=C2C%10CCCC%11C%12=CC%13=C%14CC(=C7C)C(=C3)[C@@H]%13CCC[C@H]3C7=C%13CC%15=C3C=C3C(CCCC%16C%17=CC(=C(C4)C(=C%17CC4=C%16C=C%16[C@H](CCC6)C(=C7)C(=C%13C)CC%16=C4C)C)C5CCCC1C1=CC%10=C(CC2=C8C)C(=C1C9)C)C1=CC%11=C(CC%12=C%14C)C(=C1CC3=C%15C)C)C)C";

    fn atom(element: Element) -> Atom {
        Atom::new_organic_subset(AtomSymbol::Element(element), false)
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
    fn smiles_new_and_default_create_empty_graph() {
        let smiles = Smiles::new();
        assert!(smiles.nodes().is_empty());
        assert_eq!(smiles.number_of_bonds(), 0);

        let default_smiles = Smiles::default();
        assert!(default_smiles.nodes().is_empty());
        assert_eq!(default_smiles.number_of_bonds(), 0);
    }

    #[test]
    fn bond_matrix_builder_rejects_self_loops() {
        let mut builder = BondMatrixBuilder::with_capacity(1);
        let err = builder.push_edge(0, 0, Bond::Single, None).expect_err("self-loop should fail");
        assert_eq!(err, SmilesError::SelfLoopEdge(0));
    }

    #[test]
    fn edge_key_normalizes_node_order() {
        assert_eq!(Smiles::edge_key(1, 4), (1, 4));
        assert_eq!(Smiles::edge_key(4, 1), (1, 4));
        assert_eq!(Smiles::edge_key(3, 3), (3, 3));
    }

    #[test]
    fn edge_lookup_helpers_work() {
        let ring = RingNum::try_new(1).unwrap();
        let smiles = smiles_from_edges(
            vec![atom(Element::C), atom(Element::O), atom(Element::N)],
            &[
                BondEdge::new(0, 1, Bond::Single, None),
                BondEdge::new(1, 2, Bond::Double, Some(ring)),
            ],
        );

        assert_eq!(smiles.number_of_bonds(), 2);
        assert_eq!(
            smiles.edge_for_node_pair((0, 1)),
            Some(BondEdge::new(0, 1, Bond::Single, None))
        );
        assert_eq!(
            smiles.edge_for_node_pair((1, 0)),
            Some(BondEdge::new(0, 1, Bond::Single, None))
        );
        assert_eq!(
            smiles.edge_for_node_pair((1, 2)),
            Some(BondEdge::new(1, 2, Bond::Double, Some(ring)))
        );
        assert_eq!(smiles.edge_for_node_pair((0, 2)), None);

        let edges_for_1 = smiles.edges_for_node(1);
        assert_eq!(edges_for_1.len(), 2);
        assert!(edges_for_1.contains(&BondEdge::new(1, 0, Bond::Single, None)));
        assert!(edges_for_1.contains(&BondEdge::new(1, 2, Bond::Double, Some(ring))));

        let edges_for_0 = smiles.edges_for_node(0);
        assert_eq!(edges_for_0, vec![BondEdge::new(0, 1, Bond::Single, None)]);

        let edges_for_99 = smiles.edges_for_node(99);
        assert!(edges_for_99.is_empty());
    }

    #[test]
    fn render_and_display_work_for_simple_valid_graph() {
        let smiles = smiles_from_edges(
            vec![atom(Element::C), atom(Element::O)],
            &[BondEdge::new(0, 1, Bond::Double, None)],
        );

        let rendered = smiles.render().expect("simple graph should render");
        assert_eq!(rendered, "C=O");
        assert_eq!(format!("{smiles}"), "C=O");
    }

    #[test]
    fn directional_bond_collapse_rewrites_graph_but_not_atoms() {
        let raw: Smiles = "C/C=C\\C".parse().unwrap();
        let collapsed = raw.with_directional_bonds_collapsed();

        assert_eq!(raw.nodes(), collapsed.nodes());
        assert_eq!(raw.number_of_bonds(), collapsed.number_of_bonds());
        assert_eq!(raw.edge_for_node_pair((0, 1)).unwrap().bond(), Bond::Up);
        assert_eq!(raw.edge_for_node_pair((2, 3)).unwrap().bond(), Bond::Down);
        assert_eq!(collapsed.edge_for_node_pair((0, 1)).unwrap().bond(), Bond::Single);
        assert_eq!(collapsed.edge_for_node_pair((1, 2)).unwrap().bond(), Bond::Double);
        assert_eq!(collapsed.edge_for_node_pair((2, 3)).unwrap().bond(), Bond::Single);
    }

    #[test]
    fn ring_membership_is_empty_for_acyclic_graphs() {
        let smiles: Smiles = "CCO".parse().unwrap();
        let ring_membership = smiles.try_ring_membership().unwrap();

        assert_eq!(
            ring_membership,
            RingMembership { atom_ids: Vec::new(), bond_edges: Vec::new() }
        );
        assert!(!ring_membership.contains_atom(0));
        assert!(!ring_membership.contains_edge(0, 1));
    }

    #[test]
    fn ring_membership_captures_ring_atoms_and_bonds() {
        let smiles: Smiles = "C1CCCCC1".parse().unwrap();
        let ring_membership = smiles.try_ring_membership().unwrap();

        assert_eq!(ring_membership.atom_ids(), &[0, 1, 2, 3, 4, 5]);
        assert_eq!(ring_membership.bond_edges(), &[[0, 1], [0, 5], [1, 2], [2, 3], [3, 4], [4, 5]]);
        assert!(ring_membership.contains_atom(3));
        assert!(ring_membership.contains_edge(5, 0));
        assert!(!ring_membership.contains_edge(0, 2));
    }

    #[test]
    fn symm_sssr_returns_two_fused_rings_for_naphthalene() {
        let smiles: Smiles = "c1cccc2ccccc12".parse().unwrap();

        assert_eq!(
            smiles.symm_sssr_result().cycles(),
            &[vec![0, 1, 2, 3, 4, 9], vec![4, 5, 6, 7, 8, 9]]
        );
    }

    #[test]
    fn symm_sssr_returns_cubane_faces() {
        let smiles: Smiles = "C12C3C4C1C5C2C3C45".parse().unwrap();

        assert_eq!(
            smiles.symm_sssr_result().cycles(),
            &[
                vec![0, 1, 2, 3],
                vec![0, 1, 6, 5],
                vec![0, 3, 4, 5],
                vec![1, 2, 7, 6],
                vec![2, 3, 4, 7],
                vec![4, 5, 6, 7],
            ]
        );
    }

    #[test]
    fn symm_sssr_result_reports_complete_benzene_search() {
        let smiles: Smiles = "C1=CC=CC=C1".parse().unwrap();
        let result = smiles.symm_sssr_result();

        assert!(result.status().is_complete());
        assert!(!result.status().used_fallback());
        assert!(!result.status().hit_queue_cutoff());
        assert_eq!(result.cycles(), &[vec![0, 1, 2, 3, 4, 5]]);
        assert_eq!(result.cycles(), &[vec![0, 1, 2, 3, 4, 5]]);
    }

    #[test]
    fn symm_sssr_result_exposes_fallback_for_large_polycycle_frontier_case() {
        let smiles: Smiles = LARGE_POLYCYCLE_FRONTIER_CASE.parse().unwrap();
        let result = smiles.symm_sssr_result();

        assert!(!result.status().is_complete());
        assert!(result.status().used_fallback());
        assert!(!result.status().hit_queue_cutoff());
        assert_eq!(result.cycles(), smiles.symm_sssr_result().cycles());
    }

    #[test]
    fn aromaticity_assignment_perceives_kekule_benzene() {
        let smiles: Smiles = "C1=CC=CC=C1".parse().unwrap();
        let perception = smiles.perceive_aromaticity().unwrap();
        let assignment = perception.assignment();
        let aromaticized = perception.aromaticized();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 4, 5]);
        assert_eq!(assignment.bond_edges(), &[[0, 1], [0, 5], [1, 2], [2, 3], [3, 4], [4, 5]]);
        assert!(aromaticized.nodes().iter().all(Atom::aromatic));
        assert_eq!(aromaticized.edge_for_node_pair((0, 1)).unwrap().bond(), Bond::Aromatic);
        assert_eq!(aromaticized.edge_for_node_pair((4, 5)).unwrap().bond(), Bond::Aromatic);
    }

    #[test]
    fn aromaticity_assignment_can_be_provided_by_a_custom_model() {
        struct SingleBondModel;

        impl AromaticityModel for SingleBondModel {
            fn assignment(&self, _smiles: &Smiles) -> AromaticityAssignment {
                AromaticityAssignment::new(AromaticityStatus::Complete, vec![0, 1], vec![[0, 1]])
            }
        }

        let smiles: Smiles = "CC".parse().unwrap();
        let perception = smiles.perceive_aromaticity_with(&SingleBondModel).unwrap();
        let assignment = perception.assignment();
        let aromaticized = perception.aromaticized();

        assert_eq!(assignment.atom_ids(), &[0, 1]);
        assert_eq!(assignment.bond_edges(), &[[0, 1]]);
        assert!(aromaticized.nodes().iter().all(Atom::aromatic));
        assert_eq!(aromaticized.edge_for_node_pair((0, 1)).unwrap().bond(), Bond::Aromatic);
    }

    #[test]
    fn aromaticity_policy_helpers_match_default_model() {
        let smiles: Smiles = "C1=CC=CC=C1".parse().unwrap();

        assert_eq!(
            smiles.aromaticity_assignment_for(AromaticityPolicy::RdkitDefault),
            smiles.aromaticity_assignment()
        );

        let perception = smiles.perceive_aromaticity_for(AromaticityPolicy::RdkitDefault).unwrap();
        assert_eq!(perception.status(), AromaticityStatus::Complete);
        assert_eq!(
            perception.assignment(),
            &smiles.aromaticity_assignment_for(AromaticityPolicy::RdkitDefault)
        );
        assert!(perception.aromaticized().nodes().iter().all(Atom::aromatic));

        let owned_assignment = smiles
            .perceive_aromaticity_for(AromaticityPolicy::RdkitDefault)
            .unwrap()
            .into_assignment();
        assert_eq!(owned_assignment, smiles.aromaticity_assignment());

        let owned_aromaticized = smiles
            .perceive_aromaticity_for(AromaticityPolicy::RdkitDefault)
            .unwrap()
            .into_aromaticized();
        assert!(owned_aromaticized.nodes().iter().all(Atom::aromatic));
    }

    #[test]
    fn try_with_aromaticity_assignment_accepts_valid_partial_assignment() {
        let smiles: Smiles = "CCC".parse().unwrap();
        let assignment =
            AromaticityAssignment::new(AromaticityStatus::Complete, vec![0, 1], vec![[0, 1]]);

        let aromaticized = smiles.try_with_aromaticity_assignment(&assignment).unwrap();

        assert!(aromaticized.nodes()[0].aromatic());
        assert!(aromaticized.nodes()[1].aromatic());
        assert!(!aromaticized.nodes()[2].aromatic());
        assert_eq!(aromaticized.edge_for_node_pair((0, 1)).unwrap().bond(), Bond::Aromatic);
        assert_eq!(aromaticized.edge_for_node_pair((1, 2)).unwrap().bond(), Bond::Single);
    }

    #[test]
    fn try_with_aromaticity_assignment_rejects_bonds_without_aromatic_endpoints() {
        let smiles: Smiles = "CC".parse().unwrap();
        let assignment =
            AromaticityAssignment::new(AromaticityStatus::Complete, Vec::new(), vec![[0, 1]]);

        let error = smiles.try_with_aromaticity_assignment(&assignment).unwrap_err();

        assert_eq!(
            error,
            AromaticityAssignmentApplicationError::AromaticBondMissingEndpointAtoms {
                node_a: 0,
                node_b: 1,
            }
        );
    }

    #[test]
    fn try_with_aromaticity_assignment_rejects_out_of_bounds_bond_edges() {
        let smiles: Smiles = "CC".parse().unwrap();
        let assignment =
            AromaticityAssignment::new(AromaticityStatus::Complete, vec![0, 1], vec![[0, 2]]);

        let error = smiles.try_with_aromaticity_assignment(&assignment).unwrap_err();

        assert_eq!(
            error,
            AromaticityAssignmentApplicationError::BondEdgeAtomOutOfBounds {
                node_a: 0,
                node_b: 2,
                atom_count: 2,
            }
        );
    }

    #[test]
    fn try_with_aromaticity_assignment_rejects_nonexistent_bond_edges() {
        let smiles: Smiles = "CCC".parse().unwrap();
        let assignment =
            AromaticityAssignment::new(AromaticityStatus::Complete, vec![0, 2], vec![[0, 2]]);

        let error = smiles.try_with_aromaticity_assignment(&assignment).unwrap_err();

        assert_eq!(
            error,
            AromaticityAssignmentApplicationError::MissingBondEdge { node_a: 0, node_b: 2 }
        );
    }

    #[test]
    fn try_with_aromaticity_assignment_rejects_clearing_existing_aromatic_labels() {
        let smiles: Smiles = "c1ccccc1".parse().unwrap();
        let assignment =
            AromaticityAssignment::new(AromaticityStatus::Complete, Vec::new(), Vec::new());

        let error = smiles.try_with_aromaticity_assignment(&assignment).unwrap_err();

        assert_eq!(
            error,
            AromaticityAssignmentApplicationError::WouldClearAromaticAtom { atom_id: 0 }
        );
    }

    #[test]
    fn try_with_aromaticity_assignment_rejects_out_of_bounds_atom_ids() {
        let smiles: Smiles = "CC".parse().unwrap();
        let assignment =
            AromaticityAssignment::new(AromaticityStatus::Complete, vec![2], Vec::new());

        let error = smiles.try_with_aromaticity_assignment(&assignment).unwrap_err();

        assert_eq!(
            error,
            AromaticityAssignmentApplicationError::AtomIdOutOfBounds { atom_id: 2, atom_count: 2 }
        );
    }

    #[test]
    fn try_with_aromaticity_assignment_rejects_clearing_existing_aromatic_bonds() {
        let smiles: Smiles = "c1ccccc1".parse().unwrap();
        let assignment = AromaticityAssignment::new(
            AromaticityStatus::Complete,
            vec![0, 1, 2, 3, 4, 5],
            Vec::new(),
        );

        let error = smiles.try_with_aromaticity_assignment(&assignment).unwrap_err();

        assert_eq!(
            error,
            AromaticityAssignmentApplicationError::WouldClearAromaticBond { node_a: 0, node_b: 1 }
        );
    }

    #[test]
    fn try_with_aromaticity_assignment_preserves_non_aromatic_edges_not_in_assignment() {
        let smiles: Smiles = "CCC".parse().unwrap();
        let assignment =
            AromaticityAssignment::new(AromaticityStatus::Complete, vec![0, 1, 2], vec![[0, 1]]);

        let aromaticized = smiles.try_with_aromaticity_assignment(&assignment).unwrap();

        assert_eq!(aromaticized.edge_for_node_pair((0, 1)).unwrap().bond(), Bond::Aromatic);
        assert_eq!(aromaticized.edge_for_node_pair((1, 2)).unwrap().bond(), Bond::Single);
    }

    #[test]
    fn aromaticity_assignment_preserves_custom_partial_status() {
        struct PartialModel;

        impl AromaticityModel for PartialModel {
            fn assignment(&self, _smiles: &Smiles) -> AromaticityAssignment {
                AromaticityAssignment::new(AromaticityStatus::Partial, vec![0, 1], vec![[0, 1]])
            }
        }

        let smiles: Smiles = "CC".parse().unwrap();
        let assignment = smiles.aromaticity_assignment_with(&PartialModel);

        assert_eq!(assignment.status(), AromaticityStatus::Partial);
        assert_eq!(assignment.atom_ids(), &[0, 1]);
        assert_eq!(assignment.bond_edges(), &[[0, 1]]);
    }

    #[test]
    fn perceive_aromaticity_preserves_custom_partial_status() {
        struct PartialModel;

        impl AromaticityModel for PartialModel {
            fn assignment(&self, _smiles: &Smiles) -> AromaticityAssignment {
                AromaticityAssignment::new(AromaticityStatus::Partial, vec![0, 1], vec![[0, 1]])
            }
        }

        let smiles: Smiles = "CC".parse().unwrap();
        let perception = smiles.perceive_aromaticity_with(&PartialModel).unwrap();

        assert_eq!(perception.status(), AromaticityStatus::Partial);
        assert!(perception.diagnostics().is_empty());
        assert_eq!(perception.assignment().atom_ids(), &[0, 1]);
        assert_eq!(perception.assignment().bond_edges(), &[[0, 1]]);
        assert!(perception.aromaticized().nodes().iter().all(Atom::aromatic));
        assert_eq!(
            perception.aromaticized().edge_for_node_pair((0, 1)).unwrap().bond(),
            Bond::Aromatic
        );
    }

    #[test]
    fn perceive_aromaticity_for_rdkit_default_returns_consumable_perception() {
        let smiles: Smiles = "C1=CC=CC=C1".parse().unwrap();
        let perception = smiles.perceive_aromaticity_for(AromaticityPolicy::RdkitDefault).unwrap();

        assert_eq!(perception.status(), AromaticityStatus::Complete);
        assert!(perception.diagnostics().is_empty());
        assert_eq!(perception.assignment().atom_ids(), &[0, 1, 2, 3, 4, 5]);

        let assignment = perception.into_assignment();
        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 4, 5]);

        let aromaticized = smiles
            .perceive_aromaticity_for(AromaticityPolicy::RdkitDefault)
            .unwrap()
            .into_aromaticized();
        assert!(aromaticized.nodes().iter().all(Atom::aromatic));
    }

    #[test]
    fn aromaticity_assignment_for_rdkit_default_matches_default_entrypoint() {
        let smiles: Smiles = "C1=CC=CC=C1".parse().unwrap();

        assert_eq!(
            smiles.aromaticity_assignment(),
            smiles.aromaticity_assignment_for(AromaticityPolicy::RdkitDefault)
        );
    }

    #[test]
    fn aromaticity_assignment_new_with_diagnostics_canonicalizes_inputs() {
        let assignment = AromaticityAssignment::new_with_diagnostics(
            AromaticityStatus::Partial,
            vec![3, 1, 1, 2],
            vec![[1, 3], [1, 3], [0, 2]],
            vec![
                AromaticityDiagnostic::SymmSssrHitQueueCutoff,
                AromaticityDiagnostic::SymmSssrUsedFallback,
                AromaticityDiagnostic::SymmSssrUsedFallback,
            ],
        );

        assert_eq!(assignment.atom_ids(), &[1, 2, 3]);
        assert_eq!(assignment.bond_edges(), &[[0, 2], [1, 3]]);
        assert_eq!(
            assignment.diagnostics(),
            &[
                AromaticityDiagnostic::SymmSssrUsedFallback,
                AromaticityDiagnostic::SymmSssrHitQueueCutoff,
            ]
        );
    }

    #[test]
    fn aromaticity_assignment_reports_sssr_fallback_diagnostic_for_large_polycycle_case() {
        let smiles: Smiles = LARGE_POLYCYCLE_FRONTIER_CASE.parse().unwrap();
        let assignment = smiles.aromaticity_assignment();

        assert_eq!(assignment.status(), AromaticityStatus::Partial);
        assert!(assignment.diagnostics().contains(&AromaticityDiagnostic::SymmSssrUsedFallback));
        assert!(!assignment.diagnostics().contains(&AromaticityDiagnostic::SymmSssrHitQueueCutoff));
    }

    #[test]
    fn aromaticity_assignment_perceives_kekule_pyrrole() {
        let smiles: Smiles = "C1=CNC=C1".parse().unwrap();
        let assignment = smiles.aromaticity_assignment();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 4]);
        assert_eq!(assignment.bond_edges(), &[[0, 1], [0, 4], [1, 2], [2, 3], [3, 4]]);
    }

    #[test]
    fn aromaticity_assignment_perceives_keto_pyrimidinone() {
        let smiles: Smiles = "C1=C(NC(=O)N=C1)N".parse().unwrap();
        let assignment = smiles.aromaticity_assignment();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 5, 6]);
        assert_eq!(assignment.bond_edges(), &[[0, 1], [0, 6], [1, 2], [2, 3], [3, 5], [5, 6]]);
    }

    #[test]
    fn aromaticity_assignment_handles_charged_pyridinium_systems() {
        let smiles: Smiles = "C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O)C(=O)N"
            .parse()
            .unwrap();
        let assignment = smiles.aromaticity_assignment();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 4, 5, 27, 28, 29, 30, 31, 32, 33, 34, 35]);
        assert_eq!(
            assignment.bond_edges(),
            &[
                [0, 1],
                [0, 5],
                [1, 2],
                [2, 3],
                [3, 4],
                [4, 5],
                [27, 28],
                [27, 35],
                [28, 29],
                [29, 30],
                [30, 31],
                [30, 35],
                [31, 32],
                [32, 33],
                [33, 34],
                [34, 35],
            ]
        );
    }

    #[test]
    fn aromaticity_assignment_handles_phosphinine_oxide_after_rdkit_cleanup() {
        let smiles: Smiles = "C1=CC=C2C(=C1)C=CC=P2=O".parse().unwrap();
        let assignment = smiles.aromaticity_assignment();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(
            assignment.bond_edges(),
            &[
                [0, 1],
                [0, 5],
                [1, 2],
                [2, 3],
                [3, 4],
                [3, 9],
                [4, 5],
                [4, 6],
                [6, 7],
                [7, 8],
                [8, 9]
            ]
        );
    }

    #[test]
    fn fused_symm_sssr_groups_indole_into_one_component() {
        let smiles: Smiles = "C1=CC=C2[NH]C=CC2=C1".parse().unwrap();
        let components = smiles.fused_symm_sssr_components();

        assert_eq!(components.len(), 1);
        assert_eq!(components[0].atom_ids, vec![0, 1, 2, 3, 4, 5, 6, 7, 8]);
        assert_eq!(
            components[0].bond_edges,
            vec![[0, 1], [0, 8], [1, 2], [2, 3], [3, 4], [3, 7], [4, 5], [5, 6], [6, 7], [7, 8]]
        );
    }

    #[test]
    fn aromaticity_assignment_keeps_cycloheptatriene_non_aromatic() {
        let smiles: Smiles = "C1=CC=CCC=C1".parse().unwrap();
        let assignment = smiles.aromaticity_assignment();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert!(assignment.atom_ids().is_empty());
        assert!(assignment.bond_edges().is_empty());
    }

    #[test]
    fn aromaticity_assignment_keeps_cubane_non_aromatic() {
        let smiles: Smiles = "C12C3C4C1C5C2C3C45".parse().unwrap();
        let assignment = smiles.aromaticity_assignment();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert!(assignment.atom_ids().is_empty());
        assert!(assignment.bond_edges().is_empty());
    }

    #[test]
    fn aromaticity_assignment_perceives_kekule_naphthalene() {
        let smiles: Smiles = "C1=CC=C2C=CC=CC2=C1".parse().unwrap();
        let perception = smiles.perceive_aromaticity().unwrap();
        let assignment = perception.assignment();
        let aromaticized = perception.aromaticized();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(
            assignment.bond_edges(),
            &[
                [0, 1],
                [0, 9],
                [1, 2],
                [2, 3],
                [3, 4],
                [3, 8],
                [4, 5],
                [5, 6],
                [6, 7],
                [7, 8],
                [8, 9]
            ]
        );
        assert!(aromaticized.nodes().iter().all(Atom::aromatic));
        assert_eq!(aromaticized.edge_for_node_pair((3, 8)).unwrap().bond(), Bond::Aromatic);
    }

    #[test]
    fn aromaticity_assignment_keeps_azulene_fusion_bond_non_aromatic() {
        let smiles: Smiles = "c1cccc2cccc2c1".parse().unwrap();
        let assignment = smiles.aromaticity_assignment();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(
            assignment.bond_edges(),
            &[[0, 1], [0, 9], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9]]
        );
        assert!(!assignment.contains_edge(4, 8));
    }

    #[test]
    fn aromaticity_assignment_perceives_kekule_indole() {
        let smiles: Smiles = "C1=CC=C2[NH]C=CC2=C1".parse().unwrap();
        let perception = smiles.perceive_aromaticity().unwrap();
        let assignment = perception.assignment();
        let aromaticized = perception.aromaticized();

        assert_eq!(assignment.status(), AromaticityStatus::Complete);
        assert_eq!(assignment.atom_ids(), &[0, 1, 2, 3, 4, 5, 6, 7, 8]);
        assert_eq!(
            assignment.bond_edges(),
            &[[0, 1], [0, 8], [1, 2], [2, 3], [3, 4], [3, 7], [4, 5], [5, 6], [6, 7], [7, 8]]
        );
        assert!(aromaticized.nodes().iter().all(Atom::aromatic));
        assert_eq!(aromaticized.edge_for_node_pair((3, 7)).unwrap().bond(), Bond::Aromatic);
    }

    #[test]
    fn render_smoke_from_parsed_smiles() {
        let smiles = Smiles::from_str("CC").expect("should parse");
        let rendered = smiles.render().expect("should render");
        assert_eq!(rendered, "CC");
        assert_eq!(format!("{smiles}"), "CC");
    }

    #[test]
    fn invalid_bonds_rejected() {
        assert!("-N".parse::<Smiles>().is_err());
    }

    #[test]
    fn edge_cases_branches() {
        let case_1 = "B(s)s";
        let case_1_smiles =
            case_1.parse::<Smiles>().unwrap_or_else(|_| panic!("This is a valid SMILES"));
        let case_1_smiles_rerendered = case_1_smiles.to_string();
        let case_1_smiles_reparsed = case_1_smiles_rerendered
            .parse::<Smiles>()
            .unwrap_or_else(|_| panic!("Failed to reparse B(s)s"));
        assert_eq!(case_1_smiles_reparsed.to_string(), case_1_smiles_rerendered);
    }

    #[test]
    fn branch_render_regression_non_aromatic() {
        let smiles: Smiles = "C(O)N".parse().unwrap();
        let rendered = smiles.to_string();
        assert_eq!(rendered, "C(O)N");
        let resmiles: Smiles = rendered.parse().unwrap();
        let rerendered = resmiles.to_string();
        assert_eq!(rendered, rerendered);
    }

    #[test]
    fn parse_b_s_branch_shape_is_correct() {
        let smiles: Smiles = "B(s)s".parse().unwrap();

        assert_eq!(smiles.nodes().len(), 3);
        assert_eq!(smiles.number_of_bonds(), 2);
        assert_eq!(smiles.edge_for_node_pair((0, 1)).unwrap().bond(), Bond::Single);
        assert_eq!(smiles.edge_for_node_pair((0, 2)).unwrap().bond(), Bond::Single);
        assert_eq!(smiles.edge_for_node_pair((1, 2)), None);
    }

    #[test]
    fn render_b_s_branch_should_preserve_branch_origin() {
        let smiles: Smiles = "B(s)s".parse().unwrap();
        assert_eq!(smiles.to_string(), "B(s)s");
    }

    #[test]
    fn edge_case_rings_is_invalid() {
        let err = "P1P1".parse::<Smiles>().expect_err("P1P1 should be invalid");
        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
        assert_eq!(err.start(), 3);
        assert_eq!(err.end(), 4);
    }

    #[test]
    fn edge_case_wildcard_carbon_bond() {
        let input = "*c-c";
        let smiles: Smiles = input.parse().unwrap();

        let rendered = smiles.to_string();
        let reparsed: Smiles = rendered.parse().unwrap();
        let rerendered = reparsed.to_string();

        assert_eq!(rendered, rerendered);
    }

    #[test]
    fn explicit_single_ring_closure_between_aromatic_atoms_is_preserved() {
        let smiles: Smiles = "c1-c-c-c-c-c1".parse().unwrap();
        let rendered = smiles.to_string();
        let reparsed: Smiles = rendered.parse().unwrap();
        assert_eq!(rendered, reparsed.to_string());
    }

    #[test]
    fn explicit_single_between_aromatic_atoms_is_preserved() {
        let smiles: Smiles = "*c-c".parse().unwrap();
        assert_eq!(smiles.to_string(), "*c-c");
    }

    #[test]
    fn parser_rejects_self_loop_ring() {
        let err = "C11".parse::<Smiles>().expect_err("C11 should be invalid");
        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
    }

    #[test]
    fn parser_rejects_repeated_ring_number_on_same_atom() {
        let err = "C88SC88".parse::<Smiles>().expect_err("C88SC88 should be invalid");
        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
    }

    #[test]
    fn parser_rejects_parallel_ring_bond_between_same_atoms() {
        let err = "C12CCCCC12".parse::<Smiles>().expect_err("parallel edge should be invalid");
        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
    }

    #[test]
    fn hydrogen_with_single_explicit_hydrogen_is_accepted_for_compatibility() {
        let smiles: Smiles = "[HH]".parse().expect("[HH] should be accepted for compatibility");
        assert_eq!(smiles.nodes().len(), 1);
        assert_eq!(smiles.number_of_bonds(), 0);
        assert_eq!(smiles.nodes()[0].element(), Some(Element::H));
        assert_eq!(smiles.nodes()[0].hydrogen_count(), 1);
        assert_eq!(smiles.to_string(), "[HH]");
    }

    #[test]
    fn hydrogen_with_more_than_one_explicit_hydrogen_stays_invalid() {
        let invalid =
            "[HH2]".parse::<Smiles>().expect_err("Hydrogens cannot have explicit hydrogens > 1");
        assert_eq!(invalid.smiles_error(), SmilesError::InvalidHydrogenWithExplicitHydrogensFound);
    }
}
