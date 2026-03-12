//! Module for walking the [`Smiles`] graph and tracking progress

use std::collections::HashSet;

use crate::{errors::SmilesError, smiles::Smiles, traversal::visitor_trait::Visitor};

pub fn walk<V: Visitor>(smiles: &Smiles, visitor: &mut V) -> Result<(), SmilesError> {

    Ok(())
}

