//! Wrapper module for both bracketed and unbracketed atoms

use crate::{
    atom_symbol::AtomSymbol, bracketed::bracket_atom::BracketAtom, unbracketed::UnbracketedAtom,
};

/// Enum for each variant
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Atom {
    /// [`UnbracketedAtom`] variant
    Unbracketed(UnbracketedAtom),
    /// [`BracketAtom`] variant
    Bracketed(BracketAtom),
}

impl From<UnbracketedAtom> for Atom {
    fn from(value: UnbracketedAtom) -> Self {
        Self::Unbracketed(value)
    }
}

impl From<BracketAtom> for Atom {
    fn from(value: BracketAtom) -> Self {
        Self::Bracketed(value)
    }
}

impl Atom {
    /// returns aromatic status of the atom
    #[must_use]
    pub fn aromatic(&self) -> bool {
        match self {
            Atom::Unbracketed(unbracketed_atom) => unbracketed_atom.aromatic(),
            Atom::Bracketed(bracket_atom) => bracket_atom.aromatic(),
        }
    }
    /// returns the [`AtomSymbol`]
    #[must_use]
    pub fn symbol(&self) -> AtomSymbol {
        match self {
            Atom::Unbracketed(unbracketed_atom) => unbracketed_atom.symbol(),
            Atom::Bracketed(bracket_atom) => bracket_atom.symbol(),
        }
    }
}
