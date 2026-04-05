//! Module rendering a SMILES string from the [`Smiles`] graph

use alloc::{borrow::Cow, string::String, vec, vec::Vec};
use core::{cmp::Ordering, fmt};

use crate::{
    bond::{Bond, bond_edge::BondEdge},
    errors::SmilesError,
    smiles::Smiles,
    traversal::visitor_trait::Visitor,
};

/// Structure used for implementing node visits and building output SMILES
/// `String`
pub(crate) struct RenderVisitor {
    /// Vector of output sections generated from the [`Smiles`] graph.
    sections: Vec<Cow<'static, str>>,
    /// Maps each node ID to its index in `sections`.
    node_section_index: Vec<usize>,
    /// Records where each ring label most recently ended.
    label_last_end: [usize; 100],
}

const MISSING_SECTION_INDEX: usize = usize::MAX;
const UNUSED_RING_LABEL: usize = usize::MAX;

impl RenderVisitor {
    #[must_use]
    pub(crate) fn with_capacity(number_of_nodes: usize, number_of_bonds: usize) -> Self {
        let section_capacity = number_of_nodes.saturating_mul(2).saturating_add(number_of_bonds);
        Self {
            sections: Vec::with_capacity(section_capacity),
            node_section_index: vec![MISSING_SECTION_INDEX; number_of_nodes],
            label_last_end: [UNUSED_RING_LABEL; 100],
        }
    }
    /// Builds and returns the string from the section's `String` value
    #[must_use]
    pub(crate) fn into_string(self) -> String {
        let mut output =
            String::with_capacity(self.sections.iter().map(|section| section.len()).sum());
        for section in self.sections {
            output.push_str(&section);
        }
        output
    }

    pub(crate) fn write_into_formatter(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for section in &self.sections {
            f.write_str(section.as_ref())?;
        }
        Ok(())
    }
    /// Returns the next available ring number
    ///
    /// # Errors
    /// - Returns [`SmilesError::RingNumberOverflow`] if ring number is over 99
    fn take_ring_num(&mut self, start_index: usize) -> Result<u8, SmilesError> {
        for label in 1..=99 {
            let last_end = self.label_last_end[label as usize];
            if last_end == UNUSED_RING_LABEL || last_end < start_index {
                return Ok(label);
            }
        }
        Err(SmilesError::RingNumberOverflow(100))
    }

    fn section_index_for_node(&self, node_id: usize) -> Option<usize> {
        let section_index = *self.node_section_index.get(node_id)?;
        (section_index != MISSING_SECTION_INDEX).then_some(section_index)
    }

    fn record_node_section(&mut self, node_id: usize, section_index: usize) {
        if self.node_section_index.len() <= node_id {
            self.node_section_index.resize(node_id + 1, MISSING_SECTION_INDEX);
        }
        self.node_section_index[node_id] = section_index;
    }
}

impl Visitor for RenderVisitor {
    fn enter_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError> {
        let node = smiles.nodes().get(node_id).ok_or(SmilesError::NodeIdInvalid(node_id))?;
        let index = self.sections.len();
        self.sections.push(node.rendered_cow());
        self.record_node_section(node_id, index);
        Ok(())
    }

    fn exit_node(&mut self, _smiles: &Smiles, _node_id: usize) -> Result<(), SmilesError> {
        Ok(())
    }

    fn tree_edge(&mut self, smiles: &Smiles, bond_edge: BondEdge) -> Result<(), SmilesError> {
        let symbol = if bond_edge.bond() == Bond::Single && should_render_single(smiles, bond_edge)
        {
            Bond::Single.smiles_symbol()
        } else {
            bond_edge.bond().edge_symbol()
        };
        if !symbol.is_empty() {
            self.sections.push(Cow::Borrowed(symbol));
        }
        Ok(())
    }

    fn cycle_edge(
        &mut self,
        _smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError> {
        let start_index = self.section_index_for_node(to).ok_or(SmilesError::NodeIdInvalid(to))?;
        let end_index =
            self.section_index_for_node(from).ok_or(SmilesError::NodeIdInvalid(from))?;
        let label = self.take_ring_num(start_index)?;
        match start_index.cmp(&end_index) {
            Ordering::Less => {
                let (left, right) = self.sections.split_at_mut(end_index);
                append_ring_label(left[start_index].to_mut(), label);
                append_closure_label(right[0].to_mut(), bond, label);
            }
            Ordering::Equal => {
                let section = self.sections[start_index].to_mut();
                append_ring_label(section, label);
                append_closure_label(section, bond, label);
            }
            Ordering::Greater => {
                let (left, right) = self.sections.split_at_mut(start_index);
                append_closure_label(left[end_index].to_mut(), bond, label);
                append_ring_label(right[0].to_mut(), label);
            }
        }
        self.label_last_end[label as usize] = end_index;
        Ok(())
    }
    fn open_branch(
        &mut self,
        _smiles: &Smiles,
        _from: usize,
        _to: usize,
    ) -> Result<(), SmilesError> {
        self.sections.push(Cow::Borrowed("("));
        Ok(())
    }

    fn close_branch(
        &mut self,
        _smiles: &Smiles,
        _from: usize,
        _to: usize,
    ) -> Result<(), SmilesError> {
        self.sections.push(Cow::Borrowed(")"));
        Ok(())
    }

    fn start_component(
        &mut self,
        _smiles: &Smiles,
        _root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError> {
        if component_index > 0 {
            self.sections.push(Cow::Borrowed("."));
        }
        Ok(())
    }

    fn finish_component(
        &mut self,
        _smiles: &Smiles,
        _root_id: usize,
        _component_index: usize,
    ) -> Result<(), SmilesError> {
        Ok(())
    }
}

fn should_render_single(smiles: &Smiles, bond_edge: BondEdge) -> bool {
    let (a, b) = bond_edge.vertices();
    let nodes = smiles.nodes();
    debug_assert!(a < nodes.len());
    debug_assert!(b < nodes.len());
    nodes[a].aromatic() && nodes[b].aromatic()
}

fn append_ring_label(target: &mut String, label: u8) {
    if label >= 10 {
        target.push('%');
        target.push(char::from(b'0' + (label / 10)));
        target.push(char::from(b'0' + (label % 10)));
    } else {
        target.push(char::from(b'0' + label));
    }
}

fn append_closure_label(target: &mut String, bond: Bond, label: u8) {
    target.push_str(bond.ring_closure_symbol());
    append_ring_label(target, label);
}

#[cfg(test)]
mod tests {
    use alloc::{borrow::Cow, string::ToString, vec::Vec};
    use std::str::FromStr;

    use super::{RenderVisitor, UNUSED_RING_LABEL};
    use crate::{
        bond::{Bond, bond_edge::BondEdge},
        errors::SmilesError,
        smiles::Smiles,
        traversal::visitor_trait::Visitor,
    };

    fn two_node_smiles() -> Smiles {
        "CO".parse().unwrap()
    }

    fn empty_visitor() -> RenderVisitor {
        RenderVisitor::with_capacity(0, 0)
    }

    #[test]
    fn into_string_concatenates_sections_in_order() {
        let visitor = RenderVisitor {
            sections: vec![Cow::Borrowed("C"), Cow::Borrowed("="), Cow::Borrowed("O")],
            node_section_index: Vec::new(),
            label_last_end: [UNUSED_RING_LABEL; 100],
        };

        assert_eq!(visitor.into_string(), "C=O");
    }

    #[test]
    fn enter_node_pushes_node_string_and_id() {
        let smiles = two_node_smiles();
        let mut visitor = empty_visitor();

        visitor.enter_node(&smiles, 0).unwrap();
        visitor.enter_node(&smiles, 1).unwrap();

        assert_eq!(visitor.sections.len(), 2);
        assert_eq!(visitor.sections[0], Cow::Borrowed("C"));
        assert_eq!(visitor.sections[1], Cow::Borrowed("O"));
    }

    #[test]
    fn enter_node_errors_for_invalid_node_id() {
        let smiles = two_node_smiles();
        let mut visitor = empty_visitor();

        let err = visitor.enter_node(&smiles, 99).expect_err("expected invalid node id");
        assert_eq!(err, SmilesError::NodeIdInvalid(99));
    }

    #[test]
    fn exit_node_is_ok_and_does_not_modify_sections() {
        let smiles = two_node_smiles();
        let mut visitor = empty_visitor();
        visitor.enter_node(&smiles, 0).unwrap();

        let before = visitor.sections.clone();
        visitor.exit_node(&smiles, 0).unwrap();

        assert_eq!(visitor.sections, before);
    }

    #[test]
    fn tree_edge_handles_all_bond_variants() {
        let smiles = two_node_smiles();

        let cases = [
            (Bond::Single, None),
            (Bond::Double, Some("=")),
            (Bond::Triple, Some("#")),
            (Bond::Quadruple, Some("$")),
            (Bond::Aromatic, Some(":")),
            (Bond::Up, Some("/")),
            (Bond::Down, Some("\\")),
        ];

        for (bond, expected_symbol) in cases {
            let mut visitor = empty_visitor();
            visitor.tree_edge(&smiles, BondEdge::new(0, 1, bond, None)).unwrap();

            match expected_symbol {
                None => assert!(visitor.sections.is_empty(), "single bond should add nothing"),
                Some(symbol) => {
                    assert_eq!(visitor.sections.len(), 1);
                    assert_eq!(visitor.sections[0], Cow::Borrowed(symbol));
                }
            }
        }
    }

    #[test]
    fn cycle_edge_writes_ring_labels_for_single_and_tracks_label_usage() {
        let smiles = two_node_smiles();
        let mut visitor = empty_visitor();

        visitor.enter_node(&smiles, 0).unwrap();
        visitor.enter_node(&smiles, 1).unwrap();

        visitor.cycle_edge(&smiles, 0, 1, Bond::Single).unwrap_or_else(|e| panic!("{}", e));

        assert_eq!(visitor.sections[0], Cow::Borrowed("C1"));
        assert_eq!(visitor.sections[1], Cow::Borrowed("O1"));

        assert_eq!(visitor.section_index_for_node(0), Some(0));
        assert_eq!(visitor.section_index_for_node(1), Some(1));
        assert_eq!(visitor.label_last_end[1], 0);
    }

    #[test]
    fn cycle_edge_writes_ring_labels_for_all_bond_variants() {
        let smiles = two_node_smiles();
        let mut visitor = empty_visitor();

        visitor.enter_node(&smiles, 0).unwrap();
        visitor.enter_node(&smiles, 1).unwrap();
        visitor.cycle_edge(&smiles, 0, 1, Bond::Single).unwrap();

        assert_eq!(visitor.sections[0], Cow::Borrowed("C1"));
        assert_eq!(visitor.sections[1], Cow::Borrowed("O1"));

        assert_eq!(visitor.section_index_for_node(1), Some(1));
        assert_eq!(visitor.section_index_for_node(1), Some(1));
        assert_eq!(visitor.label_last_end[1], 0);
    }

    #[test]
    fn open_and_close_branch_push_parentheses_sections() {
        let smiles = two_node_smiles();
        let mut visitor = empty_visitor();

        visitor.open_branch(&smiles, 0, 1).unwrap();
        visitor.close_branch(&smiles, 0, 1).unwrap();

        assert_eq!(visitor.sections.len(), 2);
        assert_eq!(visitor.sections[0], Cow::Borrowed("("));
        assert_eq!(visitor.sections[1], Cow::Borrowed(")"));
    }

    #[test]
    fn start_component_only_adds_dot_after_first_component() {
        let smiles = two_node_smiles();
        let mut visitor = empty_visitor();

        visitor.start_component(&smiles, 0, 0).unwrap();
        assert!(visitor.sections.is_empty());

        visitor.start_component(&smiles, 1, 1).unwrap();
        assert_eq!(visitor.sections.len(), 1);
        assert_eq!(visitor.sections[0], Cow::Borrowed("."));

        visitor.start_component(&smiles, 2, 2).unwrap();
        assert_eq!(visitor.sections.len(), 2);
        assert_eq!(visitor.sections[1], Cow::Borrowed("."));
    }

    #[test]
    fn finish_component_is_ok_and_does_not_modify_sections() {
        let smiles = two_node_smiles();
        let mut visitor = empty_visitor();
        visitor.sections.push(Cow::Borrowed("C"));

        let before = visitor.sections.clone();
        visitor.finish_component(&smiles, 0, 0).unwrap();

        assert_eq!(visitor.sections, before);
    }

    #[test]
    fn edge_case_branch_with_ring_num() {
        let input = "Ns((Ns(N)0N)N)0";

        let smiles = Smiles::from_str(input);
        assert!(smiles.is_err());
    }
    #[test]
    fn edge_case_large_oxygen_molecule() {
        let source = "C(C(F)(F)I)(CC1=CC=CC=C1N=C=NC2=CC=CC=OOF)(F)FOOC2C";
        let smiles = source.parse::<Smiles>();
        let rerender = smiles.unwrap_or_else(|e| panic!("{}", e.render(source))).to_string();
        let second_smiles = rerender.parse::<Smiles>();
        let second_rerender =
            second_smiles.unwrap_or_else(|e| panic!("{}", e.render(source))).to_string();
        assert_eq!(rerender, second_rerender);
    }
}
