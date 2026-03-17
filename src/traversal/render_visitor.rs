//! Module rendering a SMILES string from the [`Smiles`] graph

use std::collections::HashMap;

use crate::{
    bond::{Bond, bond_edge::BondEdge},
    errors::SmilesError,
    smiles::Smiles,
    traversal::visitor_trait::Visitor,
};

/// Structure used for implementing node visits and building output SMILES
/// `String`
pub struct RenderVisitor {
    /// Contains the SMILES `String` while being rendered
    output: String,
    /// [`HashMap`] used to track the ring number assigned to an edge between
    /// nodes
    ring_labels: HashMap<(usize, usize), u8>,
    /// The next available ring number label to assign for a ring
    next_ring_num: u8,
}

impl RenderVisitor {
    /// Generates a new `RenderVisitor`
    #[must_use]
    pub fn new() -> Self {
        Self { output: String::new(), ring_labels: HashMap::new(), next_ring_num: 1 }
    }
    /// Returns the built string
    #[must_use]
    pub fn into_string(self) -> String {
        self.output
    }
    /// Returns the ring label for an edge, assigning a new one if needed.
    fn ring_label_for_edge(&mut self, a: usize, b: usize) -> u8 {
        let pair = Smiles::edge_key(a, b);

        if let Some(&label) = self.ring_labels.get(&pair) {
            return label;
        }

        let label = self.next_ring_num;
        self.ring_labels.insert(pair, label);
        self.next_ring_num += 1;
        label
    }
}

impl Visitor for RenderVisitor {
    fn enter_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError> {
        if let Some(node) = smiles.node_by_id(node_id) {
            self.output.push_str(&node.to_string());
            if let Some(ring_num) = node.ring_num() {
                self.output.push_str(&ring_num.to_string());
            }
            Ok(())
        } else {
            Err(SmilesError::NodeIdInvalid(node_id))
        }
    }

    fn exit_node(&mut self, _smiles: &Smiles, _node_id: usize) -> Result<(), SmilesError> {
        Ok(())
    }

    fn tree_edge(&mut self, _smiles: &Smiles, bond_edge: BondEdge) -> Result<(), SmilesError> {
        match bond_edge.bond() {
            Bond::Single => {}
            Bond::Double => self.output.push('='),
            Bond::Triple => self.output.push('#'),
            Bond::Quadruple => self.output.push('$'),
            Bond::Aromatic => self.output.push(':'),
            Bond::Up => self.output.push('/'),
            Bond::Down => self.output.push('\\'),
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
        // let label = self.ring_label_for_edge(from, to);

        match bond {
            Bond::Single => {}
            Bond::Double => self.output.push('='),
            Bond::Triple => self.output.push('#'),
            Bond::Quadruple => self.output.push('$'),
            Bond::Aromatic => self.output.push(':'),
            Bond::Up => self.output.push('/'),
            Bond::Down => self.output.push('\\'),
        }

        // if label < 10 {
        //     self.output.push(char::from_digit(u32::from(label), 10).unwrap());
        // } else {
        //     self.output.push('%');
        //     self.output.push_str(&label.to_string());
        // }

        Ok(())
    }
    fn open_branch(
        &mut self,
        _smiles: &Smiles,
        _from: usize,
        _to: usize,
    ) -> Result<(), SmilesError> {
        self.output.push('(');
        Ok(())
    }

    fn close_branch(
        &mut self,
        _smiles: &Smiles,
        _from: usize,
        _to: usize,
    ) -> Result<(), SmilesError> {
        self.output.push(')');
        Ok(())
    }

    fn start_component(
        &mut self,
        _smiles: &Smiles,
        _root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError> {
        if component_index > 0 {
            self.output.push('.');
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
