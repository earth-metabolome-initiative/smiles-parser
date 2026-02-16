//! Module for parsing and storing the values of an atom specified in brackets -
//! `[]`.
use elements_rs::{Element, Isotope};

use crate::{
    atom_symbol::AtomSymbol,
    bracketed::{charge::Charge, chirality::Chirality, hydrogen_count::HydrogenCount},
    errors::SmilesError,
};

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Contains [`Element`] and specified meta data about an element in `[]`
pub struct BracketAtom {
    /// Bracketed elements as [`Element`]
    symbol: AtomSymbol,
    /// Parsed Isotope Mass Number Value
    isotope_mass_number: Option<u16>,
    /// If bracketed element is aromatic
    aromatic: bool,
    /// The number of Hydrogens explicitly listed or `Unspecified`
    hydrogens: HydrogenCount,
    /// The charge of the Atom, default is `0`
    charge: Charge,
    /// Esoteric potential integers presented after Element name: `[CH4:2]`.
    /// Unspecified default is 0
    class: u16,
    /// Denotes Chirality if present
    chiral: Option<Chirality>,
}

impl BracketAtom {
    /// Returns a builder for Bracket atom
    pub fn builder() -> BracketAtomBuilder {
        BracketAtomBuilder {
            bracket_atom: Self {
                symbol: AtomSymbol::default(),
                aromatic: false,
                isotope_mass_number: None,
                hydrogens: HydrogenCount::Unspecified,
                charge: Charge::default(),
                class: 0,
                chiral: None,
            },
        }
    }
    /// Returns the the [`Element`] of the bracket atom or `None` for `WildCard`
    pub fn element(&self) -> Option<Element> {
        match self.symbol {
            AtomSymbol::WildCard | AtomSymbol::Unspecified => None,
            AtomSymbol::Element(element) => Some(element),
        }
    }
    /// Returns the [`AtomSymbol`]
    pub fn symbol(&self) -> AtomSymbol {
        self.symbol
    }
    /// Returns the isotope mass number
    pub fn isotope_mass_number(&self) -> Option<u16> {
        self.isotope_mass_number
    }
    /// Returns the [`Isotope`] for the [`Element`] for the atom
    pub fn isotope(&self) -> Result<Isotope, SmilesError> {
        let element = self.element().ok_or(SmilesError::InvalidIsotope)?;
        let isotope = match self.isotope_mass_number() {
            None => element.most_abundant_isotope(),
            Some(mass) => Isotope::try_from((element, mass))?,
        };
        Ok(isotope)
    }
    /// Returns aromatic status
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }
    /// Returns the [`HydrogenCount`]
    pub fn hydrogens(&self) -> HydrogenCount {
        self.hydrogens
    }
    /// Returns the hydrogens attached or `None`
    pub fn hydrogen_count(&self) -> Option<u8> {
        match self.hydrogens {
            HydrogenCount::Unspecified => None,
            HydrogenCount::Explicit(i) => Some(i),
        }
    }
    /// Returns the [`Charge`] of the atom
    pub fn charge(&self) -> Charge {
        self.charge
    }
    /// Returns the charge as `i8`
    pub fn charge_value(&self) -> i8 {
        self.charge.get()
    }
    /// Returns the class (default is 0)
    pub fn class(&self) -> u16 {
        self.class
    }
    /// Returns the [`Chirality`] of the atom
    pub fn chiral(&self) -> Option<Chirality> {
        self.chiral
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
/// Builder structure, contains a [`BracketAtom`] that is mutable until
/// calling `build()`.
pub struct BracketAtomBuilder {
    bracket_atom: BracketAtom,
}

impl BracketAtomBuilder {
    /// Adds an isotope value
    pub fn with_isotope(mut self, iso: u16) -> Self {
        self.bracket_atom.isotope_mass_number = Some(iso);
        self
    }
    /// Adds the Atom Symbol
    pub fn with_symbol(mut self, symbol: AtomSymbol) -> Self {
        self.bracket_atom.symbol = symbol;
        self
    }
    /// Adds the aromatic value
    pub fn with_aromatic(mut self, aromatic: bool) -> Self {
        self.bracket_atom.aromatic = aromatic;
        self
    }
    /// Adds a specified [`HydrogenCount`]
    pub fn with_hydrogens(mut self, h_count: HydrogenCount) -> Self {
        self.bracket_atom.hydrogens = h_count;
        self
    }
    /// Adds a specified [`Charge`]
    pub fn with_charge(mut self, charge: Charge) -> Self {
        self.bracket_atom.charge = charge;
        self
    }
    /// Adds a designated esoteric class
    pub fn with_class(mut self, class: u16) -> Self {
        self.bracket_atom.class = class;
        self
    }
    /// Adds a specified [`Chirality`]
    pub fn with_chiral(mut self, chiral: Chirality) -> Self {
        self.bracket_atom.chiral = Some(chiral);
        self
    }
    /// Returns the [`Element`] contained in builder
    pub fn element(&self) -> Option<Element> {
        self.bracket_atom.element()
    }
    /// Returns the [`AtomSymbol`]
    pub fn symbol(&self) -> AtomSymbol {
        self.bracket_atom.symbol
    }
    /// Consumes the builder and returns the completed [`BracketAtom`]
    pub fn build(self) -> BracketAtom {
        BracketAtom {
            symbol: self.bracket_atom.symbol,
            aromatic: self.bracket_atom.aromatic,
            isotope_mass_number: self.bracket_atom.isotope_mass_number,
            hydrogens: self.bracket_atom.hydrogens,
            charge: self.bracket_atom.charge,
            class: self.bracket_atom.class,
            chiral: self.bracket_atom.chiral,
        }
    }
}
