//! Wrapper module for both bracketed and unbracketed atoms
pub mod atom_node;
pub mod atom_symbol;
pub mod bracketed;
pub mod unbracketed;

use elements_rs::Isotope;

use crate::{
    atom::{
        atom_symbol::AtomSymbol,
        bracketed::{
            BracketAtom, charge::Charge, chirality::Chirality, hydrogen_count::HydrogenCount,
        },
        unbracketed::UnbracketedAtom,
    },
    errors::SmilesError,
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
    /// returns the Chirality (if present)
    #[must_use]
    pub fn chirality(&self) -> Option<Chirality> {
        match self {
            Atom::Unbracketed(_) => None,
            Atom::Bracketed(atom) => atom.chirality(),
        }
    }
    /// returns the class (if present) of the atom
    #[must_use]
    pub fn class(&self) -> u16 {
        match self {
            Atom::Unbracketed(_) => 0,
            Atom::Bracketed(bracket_atom) => bracket_atom.class(),
        }
    }
    /// returns the [`Charge`] of the atom
    #[must_use]
    pub fn charge(&self) -> Charge {
        match self {
            Atom::Unbracketed(_) => Charge::default(),
            Atom::Bracketed(bracket_atom) => bracket_atom.charge(),
        }
    }
    /// returns the charge value as `i8`
    #[must_use]
    pub fn charge_value(&self) -> i8 {
        match self {
            Atom::Unbracketed(_) => Charge::default().get(),
            Atom::Bracketed(bracket_atom) => bracket_atom.charge_value(),
        }
    }
    /// returns the hydrogen count if present as `u8`
    #[must_use]
    pub fn hydrogen_count(&self) -> Option<u8> {
        match self {
            Atom::Unbracketed(_) => None,
            Atom::Bracketed(bracket_atom) => bracket_atom.hydrogen_count(),
        }
    }
    /// returns the [`HydrogenCount`] for the atom
    #[must_use]
    pub fn hydrogens(&self) -> HydrogenCount {
        match self {
            Atom::Unbracketed(_) => HydrogenCount::Unspecified,
            Atom::Bracketed(bracket_atom) => bracket_atom.hydrogens(),
        }
    }
    /// returns the [`Isotope`] form for the atom
    ///
    /// # Errors
    /// - Will return [`SmilesError::InvalidIsotope`] if the element is unable
    ///   to be returned
    /// - If unable to retrieve valid [`Isotope`] will return relevant
    ///   `element_rs` error
    pub fn isotope(&self) -> Result<Isotope, SmilesError> {
        match self {
            Atom::Unbracketed(unbracketed_atom) => {
                let element = unbracketed_atom.element().ok_or(SmilesError::InvalidIsotope)?;
                Ok(element.most_abundant_isotope())
            }
            Atom::Bracketed(bracket_atom) => bracket_atom.isotope(),
        }
    }
}
