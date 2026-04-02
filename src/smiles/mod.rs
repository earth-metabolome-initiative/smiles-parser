//! Represents a parsed SMILES graph.
//!
//! A [`Smiles`] value stores atoms as [`AtomNode`] values and bonds as
//! [`BondEdge`] values.
//!
//! # Examples
//!
//! SMILES strings can be parsed into a graph, inspected, and rendered back
//! into a SMILES string.
//!
//! ```rust
//! use smiles_parser::prelude::Smiles;
//!
//! let source = "CC";
//! let smiles: Smiles = source.parse()?;
//!
//! assert_eq!(smiles.nodes().len(), 2);
//! assert_eq!(smiles.edges().len(), 1);
//! assert_eq!(smiles.to_string(), "CC");
//!
//! # Ok::<(), smiles_parser::errors::SmilesErrorWithSpan>(())
//! ```
use std::fmt;

use crate::{
    atom::atom_node::AtomNode,
    bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
    errors::SmilesError,
    traversal::{render_visitor::RenderVisitor, walker::walk},
};

mod from_str;
mod implicit_hydrogens;

/// Represents a parsed SMILES graph.
#[derive(Debug, PartialEq)]
pub struct Smiles {
    atom_nodes: Vec<AtomNode>,
    bond_edges: Vec<BondEdge>,
}

impl Smiles {
    /// creates a new instance of the `Smiles` struct.
    #[must_use]
    pub fn new() -> Self {
        Self { atom_nodes: Vec::new(), bond_edges: Vec::new() }
    }
    /// Pushes an [`AtomNode`] to the `Smiles` struct.
    ///
    /// # Errors
    /// - Returns [`SmilesError::DuplicateNodeId`] if node id already exists
    pub fn push_node(&mut self, node: AtomNode) -> Result<(), SmilesError> {
        let id = node.id();
        if self.atom_nodes.get(id).is_some() {
            return Err(SmilesError::DuplicateNodeId(id));
        }
        self.atom_nodes.push(node);
        Ok(())
    }
    /// Adds a [`BondEdge`] from two nodes, includes ring number information (if
    /// present).
    ///
    /// # Errors
    /// - Returns a [`SmilesError::NodeIdInvalid`] if a node cannot be found in
    ///   the edge list
    pub fn push_edge(
        &mut self,
        node_a: usize,
        node_b: usize,
        bond: Bond,
        ring_num: Option<RingNum>,
    ) -> Result<(), SmilesError> {
        if !self.contains_node_id(node_a) {
            return Err(SmilesError::NodeIdInvalid(node_a));
        }
        if !self.contains_node_id(node_b) {
            return Err(SmilesError::NodeIdInvalid(node_b));
        }
        // reject self edges
        if node_a == node_b {
            return Err(SmilesError::SelfLoopEdge(node_a));
        }
        // reject duplicate edges
        if self.edge_for_node_pair((node_a, node_b)).is_some() {
            return Err(SmilesError::DuplicateEdge(node_a, node_b));
        }
        self.bond_edges.push(BondEdge::new(node_a, node_b, bond, ring_num));
        Ok(())
    }
    /// Returns `bool` for if the [`AtomNode`] `id` exists in the set of nodes
    /// parsed.
    fn contains_node_id(&self, id: usize) -> bool {
        self.atom_nodes.get(id).is_some()
    }
    /// Returns a slice of all [`AtomNode`] parsed in the graph.
    #[must_use]
    pub fn nodes(&self) -> &[AtomNode] {
        &self.atom_nodes
    }
    /// Returns mutable slice of [`AtomNode`] parsed in the graph.
    #[must_use]
    pub fn nodes_mut(&mut self) -> &mut [AtomNode] {
        &mut self.atom_nodes
    }
    /// Returns the node with the given `id`, if present.
    #[must_use]
    pub fn node_by_id(&self, id: usize) -> Option<&AtomNode> {
        self.atom_nodes.get(id)
    }
    /// Returns the slice of all [`BondEdge`] in the graph.
    #[must_use]
    pub fn edges(&self) -> &[BondEdge] {
        &self.bond_edges
    }
    /// Returns a normalized edge key with node IDs in ascending order. Useful
    /// for walking graph.
    ///
    /// # Parameters:
    /// - a: the first node's `id`
    /// - b: the second node's `id`
    #[must_use]
    pub fn edge_key(a: usize, b: usize) -> (usize, usize) {
        if a < b { (a, b) } else { (b, a) }
    }
    /// Returns the first [`BondEdge`] connecting the given pair of node IDs
    /// passed as a tuple.
    ///
    /// # Parameters:
    /// - nodes: A tuple of the two vertex id's.
    #[must_use]
    pub fn edge_for_node_pair(&self, nodes: (usize, usize)) -> Option<&BondEdge> {
        let target = Self::edge_key(nodes.0, nodes.1);
        self.bond_edges.iter().find(|b| {
            let (a, c) = b.vertices();
            Self::edge_key(a, c) == target
        })
    }
    /// Returns a vector of all (borrowed) [`BondEdge`] for a given [`AtomNode`]
    /// `id`.
    #[must_use]
    pub fn edges_for_node(&self, id: usize) -> Vec<&BondEdge> {
        self.bond_edges.iter().filter(|b| b.contains(id)).collect()
    }
    /// Returns mutable slice of [BondEdge].
    pub fn edges_mut(&mut self) -> &mut [BondEdge] {
        &mut self.bond_edges
    }
    /// Renders the `Smiles` graph into a valid SMILES String. Rendered `String`
    /// may differ in order and notation from input `String` but still represent
    /// the same structure.
    ///
    /// # Errors
    /// - Returns a [`SmilesError`] if the graph fails to walk
    pub fn render(&self) -> Result<String, SmilesError> {
        let mut render_visitor = RenderVisitor::new();
        walk(self, &mut render_visitor)?;
        Ok(render_visitor.into_string())
    }
}

impl Default for Smiles {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Smiles {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // need to map the smiles error to the `fmt::Error`
        let rendered_smiles = self.render().map_err(|_| fmt::Error)?;
        write!(f, "{rendered_smiles}")
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use elements_rs::Element;

    use super::Smiles;
    use crate::{
        atom::{atom_node::AtomNode, atom_symbol::AtomSymbol, unbracketed::UnbracketedAtom},
        bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
        errors::SmilesError,
    };

    fn node(id: usize, element: Element, start: usize, end: usize) -> AtomNode {
        AtomNode::new(
            UnbracketedAtom::new(AtomSymbol::Element(element), false).into(),
            id,
            start..end,
        )
    }

    #[test]
    fn smiles_new_and_default_create_empty_graph() {
        let smiles = Smiles::new();
        assert!(smiles.nodes().is_empty());
        assert!(smiles.edges().is_empty());

        let default_smiles = Smiles::default();
        assert!(default_smiles.nodes().is_empty());
        assert!(default_smiles.edges().is_empty());
    }

    #[test]
    fn push_node_adds_node_and_duplicate_id_errors() {
        let mut smiles = Smiles::new();

        let n0 = node(0, Element::C, 0, 1);
        let duplicate = node(0, Element::O, 1, 2);

        smiles.push_node(n0).expect("first node should insert");
        assert_eq!(smiles.nodes().len(), 1);
        assert_eq!(smiles.node_by_id(0), Some(&node(0, Element::C, 0, 1)));

        let err = smiles.push_node(duplicate).expect_err("duplicate id should fail");

        assert_eq!(err, SmilesError::DuplicateNodeId(0));
        assert_eq!(smiles.nodes().len(), 1);
    }

    #[test]
    fn push_edge_adds_edge_and_validates_node_ids() {
        let mut smiles = Smiles::new();
        smiles.push_node(node(0, Element::C, 0, 1)).unwrap();
        smiles.push_node(node(1, Element::O, 1, 2)).unwrap();

        smiles.push_edge(0, 1, Bond::Double, None).unwrap();

        assert_eq!(smiles.edges().len(), 1);
        assert_eq!(smiles.edges()[0], BondEdge::new(0, 1, Bond::Double, None));

        let err_a = smiles
            .push_edge(9, 1, Bond::Single, None)
            .expect_err("invalid first node id should fail");
        assert_eq!(err_a, SmilesError::NodeIdInvalid(9));

        let err_b = smiles
            .push_edge(0, 8, Bond::Single, None)
            .expect_err("invalid second node id should fail");
        assert_eq!(err_b, SmilesError::NodeIdInvalid(8));
    }

    #[test]
    fn nodes_mut_and_node_by_id_work() {
        let mut smiles = Smiles::new();
        smiles.push_node(node(0, Element::C, 0, 1)).unwrap();
        smiles.push_node(node(1, Element::O, 1, 2)).unwrap();

        assert_eq!(smiles.node_by_id(0), Some(&node(0, Element::C, 0, 1)));
        assert_eq!(smiles.node_by_id(1), Some(&node(1, Element::O, 1, 2)));
        assert_eq!(smiles.node_by_id(99), None);

        smiles.nodes_mut()[1] = node(1, Element::N, 1, 2);

        assert_eq!(smiles.node_by_id(1), Some(&node(1, Element::N, 1, 2)));
        assert_eq!(smiles.nodes()[1], node(1, Element::N, 1, 2));
    }

    #[test]
    fn edge_key_normalizes_node_order() {
        assert_eq!(Smiles::edge_key(1, 4), (1, 4));
        assert_eq!(Smiles::edge_key(4, 1), (1, 4));
        assert_eq!(Smiles::edge_key(3, 3), (3, 3));
    }

    #[test]
    fn edge_lookup_helpers_and_edges_mut_work() {
        let mut smiles = Smiles::new();
        smiles.push_node(node(0, Element::C, 0, 1)).unwrap();
        smiles.push_node(node(1, Element::O, 1, 2)).unwrap();
        smiles.push_node(node(2, Element::N, 2, 3)).unwrap();

        let ring = RingNum::try_new(1).unwrap();

        smiles.push_edge(0, 1, Bond::Single, None).unwrap();
        smiles.push_edge(1, 2, Bond::Double, Some(ring)).unwrap();

        assert_eq!(smiles.edges().len(), 2);

        assert_eq!(
            smiles.edge_for_node_pair((0, 1)),
            Some(&BondEdge::new(0, 1, Bond::Single, None))
        );
        assert_eq!(
            smiles.edge_for_node_pair((1, 0)),
            Some(&BondEdge::new(0, 1, Bond::Single, None))
        );
        assert_eq!(
            smiles.edge_for_node_pair((1, 2)),
            Some(&BondEdge::new(1, 2, Bond::Double, Some(ring)))
        );
        assert_eq!(smiles.edge_for_node_pair((0, 2)), None);

        let edges_for_1 = smiles.edges_for_node(1);
        assert_eq!(edges_for_1.len(), 2);
        assert!(edges_for_1.contains(&&BondEdge::new(0, 1, Bond::Single, None)));
        assert!(edges_for_1.contains(&&BondEdge::new(1, 2, Bond::Double, Some(ring))));

        let edges_for_0 = smiles.edges_for_node(0);
        assert_eq!(edges_for_0.len(), 1);
        assert_eq!(edges_for_0[0], &BondEdge::new(0, 1, Bond::Single, None));

        let edges_for_99 = smiles.edges_for_node(99);
        assert!(edges_for_99.is_empty());

        smiles.edges_mut()[0] = BondEdge::new(0, 1, Bond::Triple, None);
        assert_eq!(smiles.edges()[0], BondEdge::new(0, 1, Bond::Triple, None));
    }

    #[test]
    fn render_and_display_work_for_simple_valid_graph() {
        let mut smiles = Smiles::new();
        smiles.push_node(node(0, Element::C, 0, 1)).unwrap();
        smiles.push_node(node(1, Element::O, 1, 2)).unwrap();
        smiles.push_edge(0, 1, Bond::Double, None).unwrap();

        let rendered = smiles.render().expect("simple graph should render");
        assert_eq!(rendered, "C=O");

        let displayed = format!("{smiles}");
        assert_eq!(displayed, "C=O");
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
        let invalid_lead = "-N".parse::<Smiles>();
        dbg!(&invalid_lead);
        assert!(invalid_lead.is_err());
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
        assert_eq!(smiles.edges().len(), 2);

        let mut edges = smiles
            .edges()
            .iter()
            .map(|edge| {
                let (a, b) = edge.vertices();
                let key = if a < b { (a, b) } else { (b, a) };
                (key, edge.bond())
            })
            .collect::<Vec<_>>();

        edges.sort_by_key(|((a, b), _)| (*a, *b));

        assert_eq!(edges, vec![((0, 1), Bond::Single), ((0, 2), Bond::Single),]);
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
    fn hydrogen_with_explicit_hydrogen_invalid() {
        let invalid =
            "[HH]".parse::<Smiles>().expect_err("Hydrogens cannot have explicit hydrogens");
        assert_eq!(invalid.smiles_error(), SmilesError::InvalidHydrogenWithExplicitHydrogensFound);
    }
}
