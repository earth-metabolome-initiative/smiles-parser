//! Module rendering a SMILES string from the [`Smiles`] graph

use std::collections::HashMap;

use crate::{
    bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
    errors::SmilesError,
    smiles::Smiles,
    traversal::visitor_trait::Visitor,
};

/// Structure used for implementing node visits and building output SMILES
/// `String`
pub struct RenderVisitor {
    /// Vector of the outputs generated from the [`Smiles`] graph. Second value
    /// is an optional node id for sections that are nodes.
    sections: Vec<(String, Option<usize>)>,
    /// Maps each node ID to its index in `sections`.
    node_section_index: HashMap<usize, usize>,
    /// Maps each ring label to the section index where its most recent
    /// ring-closure interval ended.
    label_last_end: HashMap<u8, usize>,
}

impl RenderVisitor {
    /// Generates a new `RenderVisitor`
    #[must_use]
    pub fn new() -> Self {
        Self {
            sections: Vec::new(),
            node_section_index: HashMap::new(),
            label_last_end: HashMap::new(),
        }
    }
    /// Builds and returns the string from the section's `String` value
    #[must_use]
    pub fn into_string(self) -> String {
        let mut output = String::new();
        for (section, _) in self.sections {
            output.push_str(&section);
        }
        output
    }
    /// Returns the next available ring number
    ///
    /// # Errors
    /// - Returns [`SmilesError::RingNumberOverflow`] if ring number is over 99
    fn take_ring_num(&mut self, start_index: usize) -> Result<u8, SmilesError> {
        for label in 1..=99 {
            match self.label_last_end.get(&label) {
                Some(&last_end) if last_end >= start_index => {}
                _ => return Ok(label),
            }
        }
        Err(SmilesError::RingNumberOverflow(100))
    }
}

impl Default for RenderVisitor {
    fn default() -> Self {
        Self::new()
    }
}

impl Visitor for RenderVisitor {
    fn enter_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError> {
        if let Some(node) = smiles.node_by_id(node_id) {
            let index = self.sections.len();
            self.sections.push((node.to_string(), Some(node_id)));
            self.node_section_index.insert(node_id, index);
            Ok(())
        } else {
            Err(SmilesError::NodeIdInvalid(node_id))
        }
    }

    fn exit_node(&mut self, _smiles: &Smiles, _node_id: usize) -> Result<(), SmilesError> {
        Ok(())
    }

    fn tree_edge(&mut self, smiles: &Smiles, bond_edge: BondEdge) -> Result<(), SmilesError> {
        match bond_edge.bond() {
            Bond::Single => {
                if should_render_single(smiles, bond_edge) {
                    self.sections.push(('-'.to_string(), None));
                }
            }
            Bond::Double => self.sections.push(('='.to_string(), None)),
            Bond::Triple => self.sections.push(('#'.to_string(), None)),
            Bond::Quadruple => self.sections.push(('$'.to_string(), None)),
            Bond::Aromatic => self.sections.push((':'.to_string(), None)),
            Bond::Up => self.sections.push(('/'.to_string(), None)),
            Bond::Down => self.sections.push(('\\'.to_string(), None)),
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
        let start_index =
            *self.node_section_index.get(&to).ok_or(SmilesError::NodeIdInvalid(to))?;
        let end_index =
            *self.node_section_index.get(&from).ok_or(SmilesError::NodeIdInvalid(from))?;
        let label = self.take_ring_num(start_index)?;
        let ring_text = RingNum::try_new(label)?.to_string();
        let closure_text = match bond {
            Bond::Double => format!("={ring_text}"),
            Bond::Triple => format!("#{ring_text}"),
            Bond::Quadruple => format!("${ring_text}"),
            Bond::Single | Bond::Aromatic => ring_text.clone(),
            Bond::Up => format!("/{ring_text}"),
            Bond::Down => format!("\\{ring_text}"),
        };
        match start_index.cmp(&end_index) {
            std::cmp::Ordering::Less => {
                let (left, right) = self.sections.split_at_mut(end_index);
                left[start_index].0.push_str(&ring_text);
                right[0].0.push_str(&closure_text);
            }
            std::cmp::Ordering::Equal => {
                self.sections[start_index].0.push_str(&ring_text);
                self.sections[start_index].0.push_str(&closure_text);
            }
            std::cmp::Ordering::Greater => {
                let (left, right) = self.sections.split_at_mut(start_index);
                left[end_index].0.push_str(&closure_text);
                right[0].0.push_str(&ring_text);
            }
        }
        self.label_last_end.insert(label, end_index);
        Ok(())
    }
    fn open_branch(
        &mut self,
        _smiles: &Smiles,
        _from: usize,
        _to: usize,
    ) -> Result<(), SmilesError> {
        self.sections.push(('('.to_string(), None));
        Ok(())
    }

    fn close_branch(
        &mut self,
        _smiles: &Smiles,
        _from: usize,
        _to: usize,
    ) -> Result<(), SmilesError> {
        self.sections.push((')'.to_string(), None));
        Ok(())
    }

    fn start_component(
        &mut self,
        _smiles: &Smiles,
        _root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError> {
        if component_index > 0 {
            self.sections.push(('.'.to_string(), None));
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
    let Some(node_a) = smiles.node_by_id(a) else {
        return false;
    };
    let Some(node_b) = smiles.node_by_id(b) else {
        return false;
    };
    node_a.atom().aromatic() && node_b.atom().aromatic()
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, str::FromStr};

    use elements_rs::Element;

    use super::RenderVisitor;
    use crate::{
        atom::{atom_node::AtomNode, atom_symbol::AtomSymbol, unbracketed::UnbracketedAtom},
        bond::{Bond, bond_edge::BondEdge},
        errors::SmilesError,
        smiles::Smiles,
        traversal::visitor_trait::Visitor,
    };

    fn node(id: usize, element: Element, start: usize, end: usize) -> AtomNode {
        AtomNode::new(
            UnbracketedAtom::new(AtomSymbol::Element(element), false).into(),
            id,
            start..end,
        )
    }

    fn two_node_smiles() -> Smiles {
        let mut smiles = Smiles::new();
        smiles.push_node(node(0, Element::C, 0, 1)).unwrap();
        smiles.push_node(node(1, Element::O, 1, 2)).unwrap();
        smiles
    }

    #[test]
    fn into_string_concatenates_sections_in_order() {
        let visitor = RenderVisitor {
            sections: vec![
                ("C".to_string(), Some(0)),
                ("=".to_string(), None),
                ("O".to_string(), Some(1)),
            ],
            node_section_index: HashMap::new(),
            label_last_end: HashMap::new(),
        };

        assert_eq!(visitor.into_string(), "C=O");
    }

    #[test]
    fn enter_node_pushes_node_string_and_id() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();

        visitor.enter_node(&smiles, 0).unwrap();
        visitor.enter_node(&smiles, 1).unwrap();

        assert_eq!(visitor.sections.len(), 2);
        assert_eq!(visitor.sections[0], ("C".to_string(), Some(0)));
        assert_eq!(visitor.sections[1], ("O".to_string(), Some(1)));
    }

    #[test]
    fn enter_node_errors_for_invalid_node_id() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();

        let err = visitor.enter_node(&smiles, 99).expect_err("expected invalid node id");
        assert_eq!(err, SmilesError::NodeIdInvalid(99));
    }

    #[test]
    fn exit_node_is_ok_and_does_not_modify_sections() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();
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
            let mut visitor = RenderVisitor::new();
            visitor.tree_edge(&smiles, BondEdge::new(0, 1, bond, None)).unwrap();

            match expected_symbol {
                None => assert!(visitor.sections.is_empty(), "single bond should add nothing"),
                Some(symbol) => {
                    assert_eq!(visitor.sections.len(), 1);
                    assert_eq!(visitor.sections[0], (symbol.to_string(), None));
                }
            }
        }
    }

    #[test]
    fn cycle_edge_writes_ring_labels_for_single_and_tracks_label_usage() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();

        visitor.enter_node(&smiles, 0).unwrap();
        visitor.enter_node(&smiles, 1).unwrap();

        visitor.cycle_edge(&smiles, 0, 1, Bond::Single).unwrap_or_else(|e| panic!("{}", e));

        assert_eq!(visitor.sections[0], ("C1".to_string(), Some(0)));
        assert_eq!(visitor.sections[1], ("O1".to_string(), Some(1)));

        assert_eq!(visitor.node_section_index.get(&0), Some(&0));
        assert_eq!(visitor.node_section_index.get(&1), Some(&1));
        assert!(visitor.label_last_end.contains_key(&1));
    }

    #[test]
    fn cycle_edge_writes_ring_labels_for_all_bond_variants() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();

        visitor.enter_node(&smiles, 0).unwrap();
        visitor.enter_node(&smiles, 1).unwrap();
        visitor.cycle_edge(&smiles, 0, 1, Bond::Single).unwrap();

        assert_eq!(visitor.sections[0], ("C1".to_string(), Some(0)));
        assert_eq!(visitor.sections[1], ("O1".to_string(), Some(1)));

        assert_eq!(visitor.node_section_index.get(&1), Some(&1));
        assert_eq!(visitor.node_section_index.get(&1), Some(&1));
        assert_eq!(visitor.label_last_end.get(&1), Some(&0));
    }

    #[test]
    fn open_and_close_branch_push_parentheses_sections() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();

        visitor.open_branch(&smiles, 0, 1).unwrap();
        visitor.close_branch(&smiles, 0, 1).unwrap();

        assert_eq!(visitor.sections.len(), 2);
        assert_eq!(visitor.sections[0], ("(".to_string(), None));
        assert_eq!(visitor.sections[1], (")".to_string(), None));
    }

    #[test]
    fn start_component_only_adds_dot_after_first_component() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();

        visitor.start_component(&smiles, 0, 0).unwrap();
        assert!(visitor.sections.is_empty());

        visitor.start_component(&smiles, 1, 1).unwrap();
        assert_eq!(visitor.sections.len(), 1);
        assert_eq!(visitor.sections[0], (".".to_string(), None));

        visitor.start_component(&smiles, 2, 2).unwrap();
        assert_eq!(visitor.sections.len(), 2);
        assert_eq!(visitor.sections[1], (".".to_string(), None));
    }

    #[test]
    fn finish_component_is_ok_and_does_not_modify_sections() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();
        visitor.sections.push(("C".to_string(), Some(0)));

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
        let source = "C(C(F)(F)I)(CC1=CC=CC=C1N=C=NC2=CC=CC=OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOF)(F)FOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOC2C";
        let smiles = source.parse::<Smiles>();
        let rerender = smiles.unwrap_or_else(|e| panic!("{}", e.render(source))).to_string();
        let second_smiles = rerender.parse::<Smiles>();
        let second_rerender =
            second_smiles.unwrap_or_else(|e| panic!("{}", e.render(source))).to_string();
        assert_eq!(rerender, second_rerender);
    }
}
