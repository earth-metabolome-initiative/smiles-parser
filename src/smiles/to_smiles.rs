//! Module rendering a SMILES string from the [`Smiles`] graph

use std::collections::HashMap;

use crate::{errors::SmilesError, smiles::Smiles, traversal::visitor_trait::Visitor};

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
    /// Returns a normalized edge key with node IDs in ascending order.
    fn edge_key(a: usize, b: usize) -> (usize, usize) {
        if a < b { (a, b) } else { (b, a) }
    }
    /// Returns the ring label for an edge, assigning a new one if needed.
    fn ring_label_for_edge(&mut self, a: usize, b: usize) -> u8 {
        let pair = Self::edge_key(a, b);

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
    fn start(&mut self, smiles: &Smiles) -> Result<(), SmilesError> {
        todo!()
    }

    fn end(&mut self, smiles: &Smiles) -> Result<(), SmilesError> {
        todo!()
    }

    fn enter_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError> {
        todo!()
    }

    fn exit_node(&mut self, smiles: &Smiles, node_id: usize) -> Result<(), SmilesError> {
        todo!()
    }

    fn tree_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: crate::bond::Bond,
    ) -> Result<(), SmilesError> {
        todo!()
    }

    fn cycle_edge(
        &mut self,
        smiles: &Smiles,
        from: usize,
        to: usize,
        bond: crate::bond::Bond,
    ) -> Result<(), SmilesError> {
        todo!()
    }

    fn open_branch(&mut self, smiles: &Smiles, from: usize, to: usize) -> Result<(), SmilesError> {
        todo!()
    }

    fn close_branch(&mut self, smiles: &Smiles, from: usize, to: usize) -> Result<(), SmilesError> {
        todo!()
    }

    fn start_component(
        &mut self,
        smiles: &Smiles,
        root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError> {
        todo!()
    }

    fn finish_component(
        &mut self,
        smiles: &Smiles,
        root_id: usize,
        component_index: usize,
    ) -> Result<(), SmilesError> {
        todo!()
    }
}
