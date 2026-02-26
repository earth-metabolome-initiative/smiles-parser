//! Module for the structure of an atom as a node for use in a [`Smiles`] graph

use crate::atom::Atom;

/// Contains information about atom parsed from the SMILES string
pub struct AtomNode {
    /// Unique identifier for each node
    id: usize,
    /// Atom
    atom: Atom,
}

impl AtomNode {
    /// Creates a new node
    #[must_use]
    pub fn new(atom: Atom, id: usize) -> Self {
        Self { id, atom }
    }
    /// returns the id
    #[must_use]
    pub fn id(&self) -> usize {
        self.id
    }
    /// returns the [`Atom`]
    #[must_use]
    pub fn atom(&self) -> &Atom {
        &self.atom
    }
}