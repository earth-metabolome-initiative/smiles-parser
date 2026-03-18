//! Module rendering a SMILES string from the [`Smiles`] graph

use std::collections::HashMap;

use crate::{
    bond::{Bond, bond_edge::BondEdge, ring_num::RingNum}, errors::SmilesError, smiles::Smiles, traversal::visitor_trait::Visitor
};



/// Structure used for implementing node visits and building output SMILES
/// `String`
pub struct RenderVisitor {
    /// Vector of the outputs generated from the [`Smiles`] graph
    sections: Vec<(String, Option<usize>)>,
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
        Self { sections: Vec::new(), ring_labels: HashMap::new(), next_ring_num: 1 }
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
            self.sections.push((node.to_string(), Some(node_id)));
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
        let label = self.ring_label_for_edge(from, to);

        match bond {
            Bond::Single | Bond::Aromatic=> {}
            Bond::Double => self.sections.push(('='.to_string(), None)),
            Bond::Triple => self.sections.push(('#'.to_string(), None)),
            Bond::Quadruple => self.sections.push(('$'.to_string(), None)),
            Bond::Up => self.sections.push(('/'.to_string(), None)),
            Bond::Down => self.sections.push(('\\'.to_string(), None)),
        }

        for (output_string, id) in self.sections.iter_mut() {
            if let Some(id_val) = id {
                if *id_val == to || *id_val == from {
                    let ring_num = RingNum::try_new(label)?;
                    output_string.push_str(&ring_num.to_string());
                } 
            }
        }
        

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
