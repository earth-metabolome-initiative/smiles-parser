//! Represents a parsed SMILES graph.
//!
//! A [`Smiles`] value stores atoms as [`Atom`] values and bonds in a
//! symmetric valued sparse matrix.
//!
//! # Examples
//!
//! SMILES strings can be parsed into a graph and inspected.
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
//! # Ok::<(), smiles_parser::errors::SmilesErrorWithSpan>(())
//! ```
use alloc::{string::String, vec::Vec};
use core::fmt;

use geometric_traits::traits::{
    SizedSparseMatrix2D, SizedSparseValuedMatrixRef, SparseMatrix2D, SparseValuedMatrix2DRef,
};

use crate::{atom::Atom, bond::bond_edge::BondEdge, errors::SmilesError};

mod branches;
mod connected_components;
mod double_bond_stereo;
mod emitter;
mod from_str;
mod geometric_traits_impl;
mod implicit_hydrogens;
mod invariants;
mod neighbors;
mod refinement;
mod render_plan;
mod roots;
mod spanning_tree;
mod stereo;
mod symmetry;

pub use self::{
    connected_components::SmilesComponents,
    geometric_traits_impl::{BondEntry, BondMatrix},
};
pub(crate) use self::{geometric_traits_impl::BondMatrixBuilder, stereo::StereoNeighbor};

/// Represents a parsed SMILES graph.
#[derive(Debug)]
pub struct Smiles {
    atom_nodes: Vec<Atom>,
    bond_matrix: BondMatrix,
    parsed_stereo_neighbors: Vec<Vec<StereoNeighbor>>,
}

impl Smiles {
    /// Creates a new empty [`Smiles`] graph.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            atom_nodes: Vec::new(),
            bond_matrix: BondMatrix::default(),
            parsed_stereo_neighbors: Vec::new(),
        }
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

    /// Renders the graph back into a SMILES string.
    ///
    /// This is the main entry point for display. Rendering is split into a
    /// deterministic planning phase followed by a write-only emission phase:
    ///
    /// 1. compute atom-local invariants
    /// 2. run seeded Weisfeiler-Lehman refinement over the labeled molecular
    ///    graph
    /// 3. derive rooted symmetry tie-break classes on top of the refined
    ///    partition
    /// 4. choose one root per connected component
    /// 5. order each atom's incident edges deterministically
    /// 6. build a rooted DFS spanning forest
    /// 7. choose continuation versus branch children from subtree signatures
    /// 8. compute a final preorder for each component
    /// 9. schedule closure events and assign reusable ring labels
    /// 10. normalize stereo relative to the final traversal, with a raw
    ///     directional-single fallback for unsupported double-bond environments
    /// 11. emit a string from the finished plan
    ///
    /// The key design constraint is that ordering decisions come from
    /// graph-derived keys and explicit plan objects. The emitter does not
    /// canonicalize, search, or preview alternate strings while writing.
    ///
    /// [`fmt::Display`] uses this exact pipeline by delegating to `render()`
    /// and then writing the finished string into the formatter.
    ///
    /// # Errors
    ///
    /// This currently returns `Ok` for every in-memory graph produced by the
    /// crate and reserves `Err` for future render-time validation failures.
    pub fn render(&self) -> Result<String, SmilesError> {
        Ok(self::emitter::emit(self))
    }
}

impl Default for Smiles {
    fn default() -> Self {
        Self::new()
    }
}

impl PartialEq for Smiles {
    fn eq(&self, other: &Self) -> bool {
        self.atom_nodes == other.atom_nodes && self.bond_matrix == other.bond_matrix
    }
}

impl fmt::Display for Smiles {
    /// Formats the graph by running the full render pipeline and writing the
    /// resulting SMILES string into `f`.
    ///
    /// See [`Smiles::render`] for the full algorithm stack and the rationale
    /// behind the split between planning and emission.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rendered = self.render().map_err(|_| fmt::Error)?;
        f.write_str(&rendered)
    }
}

#[cfg(test)]
mod tests {
    use alloc::{string::ToString, vec::Vec};
    use std::str::FromStr;

    use elements_rs::Element;

    use super::{BondMatrixBuilder, Smiles};
    use crate::{
        atom::{Atom, atom_symbol::AtomSymbol},
        bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
        errors::SmilesError,
    };

    fn atom(element: Element) -> Atom {
        Atom::new_organic_subset(AtomSymbol::Element(element), false)
    }

    fn assert_render_roundtrip_stable(input: &str) {
        let smiles: Smiles = input.parse().unwrap();
        let rendered = smiles.to_string();
        let reparsed: Smiles = rendered.parse().unwrap();
        let rerendered = reparsed.to_string();
        assert_eq!(rendered, rerendered);
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
        assert_eq!(rendered, "NCO");
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
    fn render_roundtrip_stays_stable_for_non_stereo_directional_single_bond_regression() {
        assert_render_roundtrip_stable("c/o-ooNcNN-oNc");
    }

    #[test]
    fn render_roundtrip_stays_stable_for_directional_diene_regression() {
        assert_render_roundtrip_stable("C(=C/Cl)\\C=C\\Cl");
    }

    #[test]
    fn render_roundtrip_stays_stable_for_tetrahedral_ring_regression() {
        assert_render_roundtrip_stable("C[C@H]1CCCCN1N=O");
    }

    #[test]
    fn render_roundtrip_stays_stable_for_bis_tetrahedral_ring_regression() {
        assert_render_roundtrip_stable("C[C@@H]1CC[C@@H](C1)C");
    }

    #[test]
    fn render_roundtrip_stays_stable_for_cumulene_directional_regression() {
        assert_render_roundtrip_stable("C=C(=C/CO)/C#N");
    }

    #[test]
    fn render_roundtrip_stays_stable_for_polyene_boron_regression() {
        assert_render_roundtrip_stable(
            "CbbbbbbbbbbbbbbbbbbbbC/C=C\\C/C=C\\C/C=C\\CBrBrSbbbbbC#CC#CC#CC#CC#CC#CC#CC#CCCCCCCCC/C=C/C=C/CCcC/C=C/C(=C/C=C/C)CCC#CC#C",
        );
    }

    #[test]
    fn render_roundtrip_stays_stable_for_long_directional_nitroso_polycycle_regression() {
        assert_render_roundtrip_stable(
            "C=C\\1C=CC=C/OSONNNNNbcNNN:NNNNNNNNNNNC1=C-2\\C=ONNNNNbcNNN:NNNNNNNNNNNC2=C-2\\C=CCC2",
        );
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

    #[test]
    fn render_roundtrip_stays_stable_for_pubchem_high_ring_label_regression() {
        assert_render_roundtrip_stable(concat!(
            "C1C2C3C4C5C6C7C8CC9C88C77C66C55C44C33C22C1C1C22C33C44C55C66C77C88C9C9C88",
            "C77C66C55C44C33C22C1C1C22C33C44C55C66C77C88C9C9C88C77C66C55C44C33C22C1C1",
            "C22C33C44C55C66C77C88C9C9C88C77C66C55C44C33C22C1C1C22C33C44C55C66C77C88C",
            "9C9C88C77C66C55C44C33C22C1C1C22C33C44C55C66C77C88C9C9C88C77C66C55C44C33C",
            "22C1C1C22C33C44C55C66C77C88C9C9C88C77C66C55C44C33C22C1C1C22C33C44C55C66C",
            "77C88C9C9C88C77C66C55C44C33C22C1C1C22C33C44C55C66C77C88C9C9C88C77C66C55C",
            "44C33C22C1C1C22C33C44C55C66C77C88C9C9C88C77C66C55C44C33C22C1C1C22C33C44C",
            "55C66C77C88C9C9C88C77C66C55C44C33C22C1C1C22C33C44C55C66C77C88C9C9C88C77C",
            "66C55C44C33C22C1C1C22C33C44C55C66C77C88C9C9C88C77C66C55C44C33C22C1C1C22C",
            "33C44C55C66C77C88C9C9C88C77C66C55C44C33C22C1C1C22C33C44C55C66C77C88C9C9C",
            "88C77C66C55C44C33C22C1C1C22C33C44C55C66C77C88C9C9C88C77C66C55C44C33C22C1",
            "CC2C3C4C5C6C7C8C9",
        ));
    }

    #[test]
    fn render_roundtrip_stays_stable_for_quaternary_cage_hydroxide_regression() {
        assert_render_roundtrip_stable(
            "C1N2CN3CN1C[N+](C2)(C3)CO.C1N2CN3CN1C[N+](C2)(C3)CO.O.O.O.O.[OH-].[OH-].[OH-].[OH-]",
        );
    }
}
