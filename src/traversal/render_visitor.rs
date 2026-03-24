//! Module rendering a SMILES string from the [`Smiles`] graph

use std::collections::BTreeSet;

use crate::{
    bond::{Bond, bond_edge::BondEdge, ring_num::RingNum},
    errors::SmilesError,
    smiles::Smiles,
    traversal::visitor_trait::Visitor,
};

/// Structure used for implementing node visits and building output SMILES
/// `String`
pub struct RenderVisitor {
    /// Vector of the outputs generated from the [`Smiles`] graph
    sections: Vec<(String, Option<usize>)>,
    /// A pool of ring numbers that can be used
    available_ring_nums: BTreeSet<u8>,
}

impl RenderVisitor {
    /// Generates a new `RenderVisitor`
    #[must_use]
    pub fn new() -> Self {
        Self { sections: Vec::new(), available_ring_nums: (0..=99).collect() }
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
    fn take_ring_num(&mut self) -> Result<u8, SmilesError> {
        let Some(&ring_num) = self.available_ring_nums.iter().next() else {
            return Err(SmilesError::RingNumberOverflow(100));
        };
        self.available_ring_nums.remove(&ring_num);
        Ok(ring_num)
    }
    /// Recycles a ring number after its finished being used
    fn release_ring_num(&mut self, ring_num: u8) {
        self.available_ring_nums.insert(ring_num);
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
            self.sections.push((node.to_string(), Some(node_id)));
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
                    self.sections.push(('-'.to_string(), None))
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
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: Bond,
    ) -> Result<(), SmilesError> {
        let label = self.take_ring_num()?;
        let ring_text = RingNum::try_new(label)?.to_string();
        let closure_text = match bond {
            Bond::Single => {
                let Some(node_from) = smiles.node_by_id(from) else {
                    return Err(SmilesError::NodeIdInvalid(from));
                };
                let Some(node_to) = smiles.node_by_id(to) else {
                    return Err(SmilesError::NodeIdInvalid(to));
                };
                if node_from.atom().aromatic() && node_to.atom().aromatic() {
                    format!("-{ring_text}")
                } else {
                    ring_text.clone()
                }
            }
            Bond::Double => format!("={ring_text}"),
            Bond::Triple => format!("#{ring_text}"),
            Bond::Quadruple => format!("${ring_text}"),
            Bond::Aromatic => ring_text.clone(),
            Bond::Up => format!("/{ring_text}"),
            Bond::Down => format!("\\{ring_text}"),
        };

        for (output_string, id) in &mut self.sections {
            if let Some(id_val) = id {
                if *id_val == to {
                    output_string.push_str(&ring_text);
                } else if *id_val == from {
                    output_string.push_str(&closure_text);
                }
            }
        }
        self.release_ring_num(label);
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
    use std::collections::BTreeSet;

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
    fn new_and_default_initialize_empty_sections_and_full_ring_pool() {
        let visitor = RenderVisitor::new();
        assert!(visitor.sections.is_empty());
        assert_eq!(visitor.available_ring_nums.len(), 100);
        assert!(visitor.available_ring_nums.contains(&0));
        assert!(visitor.available_ring_nums.contains(&99));

        let default_visitor = RenderVisitor::default();
        assert!(default_visitor.sections.is_empty());
        assert_eq!(default_visitor.available_ring_nums.len(), 100);
    }

    #[test]
    fn into_string_concatenates_sections_in_order() {
        let visitor = RenderVisitor {
            sections: vec![
                ("C".to_string(), Some(0)),
                ("=".to_string(), None),
                ("O".to_string(), Some(1)),
            ],
            available_ring_nums: (0..=99).collect(),
        };

        assert_eq!(visitor.into_string(), "C=O");
    }

    #[test]
    fn take_ring_num_returns_lowest_available_and_release_ring_num_restores_it() {
        let mut visitor = RenderVisitor::new();

        let first = visitor.take_ring_num().unwrap();
        assert_eq!(first, 0);
        assert!(!visitor.available_ring_nums.contains(&0));
        assert_eq!(visitor.available_ring_nums.len(), 99);

        let second = visitor.take_ring_num().unwrap();
        assert_eq!(second, 1);
        assert!(!visitor.available_ring_nums.contains(&1));
        assert_eq!(visitor.available_ring_nums.len(), 98);

        visitor.release_ring_num(first);
        assert!(visitor.available_ring_nums.contains(&0));
        assert_eq!(visitor.available_ring_nums.len(), 99);

        let recycled = visitor.take_ring_num().unwrap();
        assert_eq!(recycled, 0);
    }

    #[test]
    fn take_ring_num_errors_when_pool_is_empty() {
        let mut visitor =
            RenderVisitor { sections: Vec::new(), available_ring_nums: BTreeSet::default() };

        let err = visitor.take_ring_num().expect_err("expected ring number overflow");
        assert_eq!(err, SmilesError::RingNumberOverflow(100));
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
    fn cycle_edge_writes_ring_labels_for_single_and_releases_number() {
        let smiles = two_node_smiles();
        let mut visitor = RenderVisitor::new();

        visitor.enter_node(&smiles, 0).unwrap();
        visitor.enter_node(&smiles, 1).unwrap();

        visitor.cycle_edge(&smiles, 0, 1, Bond::Single).unwrap();

        assert_eq!(visitor.sections[0], ("C0".to_string(), Some(0)));
        assert_eq!(visitor.sections[1], ("O0".to_string(), Some(1)));
        assert!(visitor.available_ring_nums.contains(&0));
        assert_eq!(visitor.available_ring_nums.len(), 100);
    }

    #[test]
    fn cycle_edge_writes_ring_labels_for_all_bond_variants() {
        let smiles = two_node_smiles();

        let cases = [
            (Bond::Single, "0", "0"),
            (Bond::Aromatic, "0", "0"),
            (Bond::Double, "=0", "0"),
            (Bond::Triple, "#0", "0"),
            (Bond::Quadruple, "$0", "0"),
            (Bond::Up, "/0", "0"),
            (Bond::Down, "\\0", "0"),
        ];

        for (bond, from_suffix, to_suffix) in cases {
            let mut visitor = RenderVisitor::new();
            visitor.enter_node(&smiles, 0).unwrap();
            visitor.enter_node(&smiles, 1).unwrap();

            visitor.cycle_edge(&smiles, 0, 1, bond).unwrap();

            assert_eq!(visitor.sections[0], (format!("C{from_suffix}"), Some(0)));
            assert_eq!(visitor.sections[1], (format!("O{to_suffix}"), Some(1)));
            assert!(visitor.available_ring_nums.contains(&0));
        }
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
}
