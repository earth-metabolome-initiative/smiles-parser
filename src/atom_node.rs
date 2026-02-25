//! Module for the structure of an atom as a node for use in a [`Smiles`] graph

use crate::{atom::Atom, atom_symbol::AtomSymbol};

/// Contains information about atom parsed from the SMILES string
pub struct AtomNode {
    /// Unique identifier for each node
    id: usize,
    /// Atom
    atom: Atom,
}
