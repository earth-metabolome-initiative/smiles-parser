//! Aromaticity policy and assignment utilities for [`Smiles`](super::Smiles).

use alloc::vec::Vec;

use geometric_traits::traits::SparseValuedMatrixRef;
use thiserror::Error;

use super::Smiles;
use crate::bond::Bond;

mod rdkit_default;

/// Aromaticity-model interface for [`Smiles`](super::Smiles).
///
/// Implementations take a parsed graph and produce a canonicalized aromaticity
/// assignment. The assignment can then be inspected directly or applied back to
/// a transformed [`Smiles`](super::Smiles) graph.
pub trait AromaticityModel {
    /// Perceives aromatic atoms and bonds for the provided graph.
    fn assignment(&self, smiles: &Smiles) -> AromaticityAssignment;
}

/// Named aromaticity-policy presets exposed by the crate.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum AromaticityPolicy {
    /// The current `RDKit` default aromaticity model target.
    RdkitDefault,
}

impl AromaticityModel for AromaticityPolicy {
    fn assignment(&self, smiles: &Smiles) -> AromaticityAssignment {
        match self {
            Self::RdkitDefault => RdkitDefaultAromaticity.assignment(smiles),
        }
    }
}

/// The `RDKit` default aromaticity-model implementation.
#[derive(Debug, Copy, Clone, Default, PartialEq, Eq, Hash)]
pub struct RdkitDefaultAromaticity;

/// Completeness of an aromaticity assignment.
#[non_exhaustive]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum AromaticityStatus {
    /// The assignment completed without encountering unsupported cases.
    Complete,
    /// The assignment completed, but at least one internal stage used a
    /// fallback, hit a search cutoff, or skipped unsupported work.
    Partial,
    /// The assignment could not be completed for the requested model.
    Unsupported,
}

/// Ring-family category used by aromaticity diagnostics.
#[non_exhaustive]
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticityRingFamilyKind {
    /// A single simple cycle.
    SimpleCycle,
    /// A fused family composed of multiple cycles.
    FusedComponent,
}

/// Diagnostic reason explaining why an aromaticity assignment is partial or
/// unsupported.
#[non_exhaustive]
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticityDiagnostic {
    /// Symmetrized SSSR had to fall back to approximate ring finding.
    SymmSssrUsedFallback,
    /// Symmetrized SSSR hit the explicit BFS work budget.
    SymmSssrHitQueueCutoff,
    /// A fused-family search stopped at the configured exact-search budget.
    FusedSubsystemBudgetExceeded {
        /// Number of rings in the fused family.
        ring_count: usize,
        /// Maximum fused-family ring count searched exactly.
        max_fused_subsystem_rings: usize,
        /// Maximum two-ring frontier size searched exactly.
        max_ring_combination_search: usize,
    },
    /// A ring family was skipped as unsupported by the current model.
    UnsupportedRingFamily {
        /// Ring-family category that was skipped.
        kind: AromaticityRingFamilyKind,
        /// Number of member cycles in the skipped family.
        ring_count: usize,
        /// Number of unique atoms in the skipped family.
        atom_count: usize,
        /// Number of unique bonds in the skipped family.
        bond_count: usize,
    },
}

/// Error raised while applying an aromaticity assignment back onto a
/// [`Smiles`](super::Smiles) graph.
///
/// Applying an assignment is only lossless when the source graph does not
/// require clearing pre-existing aromatic labels. Once a graph already contains
/// [`Bond::Aromatic`] edges, the original non-aromatic bond orders are no
/// longer recoverable from the transformed graph alone.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Error)]
pub enum AromaticityAssignmentApplicationError {
    /// The assignment references an atom id outside the graph.
    #[error("aromaticity assignment atom index out of bounds: {atom_id} (atom count {atom_count})")]
    AtomIdOutOfBounds {
        /// Referenced atom id.
        atom_id: usize,
        /// Number of atoms in the graph.
        atom_count: usize,
    },
    /// The assignment references a bond edge whose endpoints are outside the
    /// graph.
    #[error(
        "aromaticity assignment bond edge references invalid atom index: [{node_a}, {node_b}] (atom count {atom_count})"
    )]
    BondEdgeAtomOutOfBounds {
        /// First endpoint of the referenced bond edge.
        node_a: usize,
        /// Second endpoint of the referenced bond edge.
        node_b: usize,
        /// Number of atoms in the graph.
        atom_count: usize,
    },
    /// The assignment references a bond edge that does not exist in the graph.
    #[error("aromaticity assignment bond edge does not exist in graph: [{node_a}, {node_b}]")]
    MissingBondEdge {
        /// First endpoint of the missing bond edge.
        node_a: usize,
        /// Second endpoint of the missing bond edge.
        node_b: usize,
    },
    /// Aromatic bond assignments require both endpoint atoms to be aromatic.
    #[error(
        "aromaticity assignment bond edge [{node_a}, {node_b}] requires both endpoint atoms to be aromatic"
    )]
    AromaticBondMissingEndpointAtoms {
        /// First endpoint of the inconsistent aromatic bond edge.
        node_a: usize,
        /// Second endpoint of the inconsistent aromatic bond edge.
        node_b: usize,
    },
    /// Applying the assignment would require clearing an already-aromatic atom.
    #[error(
        "cannot apply aromaticity assignment because it would clear pre-existing aromatic atom {atom_id}"
    )]
    WouldClearAromaticAtom {
        /// Already-aromatic atom that the assignment excludes.
        atom_id: usize,
    },
    /// Applying the assignment would require clearing an already-aromatic bond.
    #[error(
        "cannot apply aromaticity assignment because it would clear pre-existing aromatic bond [{node_a}, {node_b}]"
    )]
    WouldClearAromaticBond {
        /// First endpoint of the already-aromatic bond that the assignment
        /// excludes.
        node_a: usize,
        /// Second endpoint of the already-aromatic bond that the assignment
        /// excludes.
        node_b: usize,
    },
}

/// Canonicalized aromatic atom and bond assignment for a
/// [`Smiles`](super::Smiles) graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AromaticityAssignment {
    status: AromaticityStatus,
    atom_ids: Vec<usize>,
    bond_edges: Vec<[usize; 2]>,
    diagnostics: Vec<AromaticityDiagnostic>,
}

impl AromaticityAssignment {
    /// Creates a canonicalized aromaticity assignment.
    #[must_use]
    pub fn new(
        status: AromaticityStatus,
        atom_ids: Vec<usize>,
        bond_edges: Vec<[usize; 2]>,
    ) -> Self {
        Self::new_with_diagnostics(status, atom_ids, bond_edges, Vec::new())
    }

    /// Creates a canonicalized aromaticity assignment with diagnostic reasons.
    #[must_use]
    pub fn new_with_diagnostics(
        status: AromaticityStatus,
        atom_ids: Vec<usize>,
        bond_edges: Vec<[usize; 2]>,
        mut diagnostics: Vec<AromaticityDiagnostic>,
    ) -> Self {
        let mut atom_ids = atom_ids;
        let mut bond_edges = bond_edges;
        atom_ids.sort_unstable();
        atom_ids.dedup();
        bond_edges.sort_unstable();
        bond_edges.dedup();
        diagnostics.sort_unstable();
        diagnostics.dedup();
        Self { status, atom_ids, bond_edges, diagnostics }
    }

    /// Returns the aromaticity assignment status.
    #[inline]
    #[must_use]
    pub fn status(&self) -> AromaticityStatus {
        self.status
    }

    /// Returns the aromatic atom ids.
    #[inline]
    #[must_use]
    pub fn atom_ids(&self) -> &[usize] {
        &self.atom_ids
    }

    /// Returns the aromatic bond edges.
    #[inline]
    #[must_use]
    pub fn bond_edges(&self) -> &[[usize; 2]] {
        &self.bond_edges
    }

    /// Returns the diagnostic reasons attached to the assignment.
    #[inline]
    #[must_use]
    pub fn diagnostics(&self) -> &[AromaticityDiagnostic] {
        &self.diagnostics
    }

    /// Returns whether the given atom id is assigned aromatic.
    #[inline]
    #[must_use]
    pub fn contains_atom(&self, atom_id: usize) -> bool {
        self.atom_ids.binary_search(&atom_id).is_ok()
    }

    /// Returns whether the given edge is assigned aromatic.
    #[inline]
    #[must_use]
    pub fn contains_edge(&self, node_a: usize, node_b: usize) -> bool {
        let edge = if node_a < node_b { [node_a, node_b] } else { [node_b, node_a] };
        self.bond_edges.binary_search(&edge).is_ok()
    }

    /// Validates that this assignment can be losslessly applied to the
    /// provided graph.
    ///
    /// This rejects inconsistent assignments and assignments that would require
    /// clearing pre-existing aromatic labels from the source graph.
    ///
    /// # Errors
    /// Returns an [`AromaticityAssignmentApplicationError`] if the assignment
    /// is inconsistent with the graph or would require clearing information the
    /// graph no longer stores.
    pub fn validate_for(
        &self,
        smiles: &Smiles,
    ) -> Result<(), AromaticityAssignmentApplicationError> {
        let atom_count = smiles.nodes().len();

        for &atom_id in &self.atom_ids {
            if atom_id >= atom_count {
                return Err(AromaticityAssignmentApplicationError::AtomIdOutOfBounds {
                    atom_id,
                    atom_count,
                });
            }
        }

        for &[node_a, node_b] in &self.bond_edges {
            if node_a >= atom_count || node_b >= atom_count {
                return Err(AromaticityAssignmentApplicationError::BondEdgeAtomOutOfBounds {
                    node_a,
                    node_b,
                    atom_count,
                });
            }
            if smiles.edge_for_node_pair((node_a, node_b)).is_none() {
                return Err(AromaticityAssignmentApplicationError::MissingBondEdge {
                    node_a,
                    node_b,
                });
            }
            if !self.contains_atom(node_a) || !self.contains_atom(node_b) {
                return Err(
                    AromaticityAssignmentApplicationError::AromaticBondMissingEndpointAtoms {
                        node_a,
                        node_b,
                    },
                );
            }
        }

        for (atom_id, atom) in smiles.nodes().iter().enumerate() {
            if atom.aromatic() && !self.contains_atom(atom_id) {
                return Err(AromaticityAssignmentApplicationError::WouldClearAromaticAtom {
                    atom_id,
                });
            }
        }

        for ((row, column), entry) in smiles.bond_matrix().sparse_entries() {
            if row < column && entry.bond() == Bond::Aromatic && !self.contains_edge(row, column) {
                return Err(AromaticityAssignmentApplicationError::WouldClearAromaticBond {
                    node_a: row,
                    node_b: column,
                });
            }
        }

        Ok(())
    }
}

/// A completed aromaticity perception pass together with the transformed
/// graph it produces.
#[derive(Debug, PartialEq)]
pub struct AromaticityPerception {
    assignment: AromaticityAssignment,
    aromaticized: Smiles,
}

impl AromaticityPerception {
    /// Returns the aromaticity assignment.
    #[inline]
    #[must_use]
    pub fn assignment(&self) -> &AromaticityAssignment {
        &self.assignment
    }

    /// Returns the aromaticity status.
    #[inline]
    #[must_use]
    pub fn status(&self) -> AromaticityStatus {
        self.assignment.status()
    }

    /// Returns the diagnostic reasons attached to the assignment.
    #[inline]
    #[must_use]
    pub fn diagnostics(&self) -> &[AromaticityDiagnostic] {
        self.assignment.diagnostics()
    }

    /// Returns the transformed graph with aromatic labels applied.
    #[inline]
    #[must_use]
    pub fn aromaticized(&self) -> &Smiles {
        &self.aromaticized
    }

    /// Consumes the perception and returns the aromaticity assignment.
    #[inline]
    #[must_use]
    pub fn into_assignment(self) -> AromaticityAssignment {
        self.assignment
    }

    /// Consumes the perception and returns the aromaticized graph.
    #[inline]
    #[must_use]
    pub fn into_aromaticized(self) -> Smiles {
        self.aromaticized
    }
}

impl Smiles {
    /// Returns the current default aromaticity assignment.
    ///
    /// This delegates to [`RdkitDefaultAromaticity`].
    ///
    /// The implementation targets `RDKit`'s default aromaticity model and is
    /// validated against the frozen fixture corpora plus the extracted whole-
    /// `PubChem` `RDKit`-default aromatic corpus used by the test suite.
    #[inline]
    #[must_use]
    pub fn aromaticity_assignment(&self) -> AromaticityAssignment {
        self.aromaticity_assignment_with(&RdkitDefaultAromaticity)
    }

    /// Returns the aromaticity assignment produced by the provided named
    /// policy preset.
    #[inline]
    #[must_use]
    pub fn aromaticity_assignment_for(&self, policy: AromaticityPolicy) -> AromaticityAssignment {
        self.aromaticity_assignment_with(&policy)
    }

    /// Returns the aromaticity assignment produced by the provided model.
    #[inline]
    #[must_use]
    pub fn aromaticity_assignment_with<M: AromaticityModel>(
        &self,
        model: &M,
    ) -> AromaticityAssignment {
        model.assignment(self)
    }

    /// Returns the current default aromaticity perception result.
    ///
    /// The returned value keeps both the aromaticity assignment and the
    /// transformed graph together so callers can inspect completeness status
    /// before discarding it.
    ///
    /// # Errors
    /// Returns an [`AromaticityAssignmentApplicationError`] if the assignment
    /// is inconsistent with the graph or would require clearing pre-existing
    /// aromatic labels from the source graph.
    pub fn perceive_aromaticity(
        &self,
    ) -> Result<AromaticityPerception, AromaticityAssignmentApplicationError> {
        self.perceive_aromaticity_with(&RdkitDefaultAromaticity)
    }

    /// Returns the aromaticity perception result for the provided named policy
    /// preset.
    ///
    /// # Errors
    /// Returns an [`AromaticityAssignmentApplicationError`] if the assignment
    /// is inconsistent with the graph or would require clearing pre-existing
    /// aromatic labels from the source graph.
    pub fn perceive_aromaticity_for(
        &self,
        policy: AromaticityPolicy,
    ) -> Result<AromaticityPerception, AromaticityAssignmentApplicationError> {
        self.perceive_aromaticity_with(&policy)
    }

    /// Returns the aromaticity perception result for the provided model.
    ///
    /// # Errors
    /// Returns an [`AromaticityAssignmentApplicationError`] if the assignment
    /// is inconsistent with the graph or would require clearing pre-existing
    /// aromatic labels from the source graph.
    pub fn perceive_aromaticity_with<M: AromaticityModel>(
        &self,
        model: &M,
    ) -> Result<AromaticityPerception, AromaticityAssignmentApplicationError> {
        let assignment = self.aromaticity_assignment_with(model);
        let aromaticized = self.try_with_aromaticity_assignment(&assignment)?;
        Ok(AromaticityPerception { assignment, aromaticized })
    }

    /// Returns a new graph with the provided aromaticity assignment applied to
    /// atom and bond labels.
    ///
    /// # Errors
    /// Returns an [`AromaticityAssignmentApplicationError`] if the assignment
    /// is inconsistent with the graph or would require clearing pre-existing
    /// aromatic labels from the source graph.
    pub fn try_with_aromaticity_assignment(
        &self,
        assignment: &AromaticityAssignment,
    ) -> Result<Self, AromaticityAssignmentApplicationError> {
        assignment.validate_for(self)?;

        let atom_nodes =
            self.atom_nodes
                .iter()
                .copied()
                .enumerate()
                .map(|(atom_id, atom)| {
                    if assignment.contains_atom(atom_id) { atom.with_aromatic(true) } else { atom }
                })
                .collect::<Vec<_>>();

        let bond_matrix = super::BondMatrix::from_sorted_upper_triangular_entries(
            atom_nodes.len(),
            self.bond_matrix.sparse_entries().filter_map(|((row, column), entry)| {
                (row < column).then_some((
                    row,
                    column,
                    if assignment.contains_edge(row, column) {
                        entry.with_bond(Bond::Aromatic)
                    } else {
                        *entry
                    },
                ))
            }),
        )
        .unwrap_or_else(|_| unreachable!("existing bond matrix entries are already valid"));

        Ok(Self { atom_nodes, bond_matrix })
    }
}
