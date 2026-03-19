//! Module rendering a SMILES string from the [`Smiles`] graph

use std::collections::BTreeSet;

use crate::{
    bond::{
        Bond,
        bond_edge::BondEdge,
        ring_num::RingNum,
    },
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
        let label = self.take_ring_num()?;
        let ring_text = RingNum::try_new(label)?.to_string();
        let closure_text = match bond {
            Bond::Single | Bond::Aromatic => ring_text.clone(),
            Bond::Double => format!("={ring_text}"),
            Bond::Triple => format!("#{ring_text}"),
            Bond::Quadruple => format!("${ring_text}"),
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
