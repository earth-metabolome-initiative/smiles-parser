//! Module for parsing and storing the values of an atom specified in brackets -
//! `[]`.

pub mod charge;
pub mod chirality;
pub mod hydrogen_count;
use std::fmt;

use elements_rs::{Element, Isotope};

use crate::{
    atom::{
        atom_symbol::AtomSymbol,
        bracketed::{charge::Charge, chirality::Chirality, hydrogen_count::HydrogenCount},
    },
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
    chirality: Option<Chirality>,
}

impl BracketAtom {
    /// Returns a builder for Bracket atom
    #[must_use]
    pub fn builder() -> BracketAtomBuilder {
        BracketAtomBuilder {
            bracket_atom: Self {
                symbol: AtomSymbol::default(),
                aromatic: false,
                isotope_mass_number: None,
                hydrogens: HydrogenCount::Unspecified,
                charge: Charge::default(),
                class: 0,
                chirality: None,
            },
        }
    }
    /// Returns the the [`Element`] of the bracket atom or `None` for `WildCard`
    #[must_use]
    pub fn element(&self) -> Option<Element> {
        match self.symbol {
            AtomSymbol::WildCard | AtomSymbol::Unspecified => None,
            AtomSymbol::Element(element) => Some(element),
        }
    }
    /// Returns the [`AtomSymbol`]
    #[must_use]
    pub fn symbol(&self) -> AtomSymbol {
        self.symbol
    }
    /// Returns the isotope mass number
    #[must_use]
    pub fn isotope_mass_number(&self) -> Option<u16> {
        self.isotope_mass_number
    }
    /// Returns the [`Isotope`] for the [`Element`] for the atom
    ///
    /// # Errors
    /// - Returns [`SmilesError::InvalidIsotope`] the element fails to parse
    /// - Returns `element_rs` error if unable to return correct isotope
    pub fn isotope(&self) -> Result<Isotope, SmilesError> {
        let element = self.element().ok_or(SmilesError::InvalidIsotope)?;
        let isotope = match self.isotope_mass_number() {
            None => element.most_abundant_isotope(),
            Some(mass) => Isotope::try_from((element, mass))?,
        };
        Ok(isotope)
    }
    /// Returns aromatic status
    #[must_use]
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }
    /// Returns the [`HydrogenCount`]
    #[must_use]
    pub fn hydrogens(&self) -> HydrogenCount {
        self.hydrogens
    }
    /// Returns the hydrogens attached or `None`
    #[must_use]
    pub fn hydrogen_count(&self) -> Option<u8> {
        match self.hydrogens {
            HydrogenCount::Unspecified => None,
            HydrogenCount::Explicit(i) => Some(i),
        }
    }
    /// Returns the [`Charge`] of the atom
    #[must_use]
    pub fn charge(&self) -> Charge {
        self.charge
    }
    /// Returns the charge as `i8`
    #[must_use]
    pub fn charge_value(&self) -> i8 {
        self.charge.get()
    }
    /// Returns the class (default is 0)
    #[must_use]
    pub fn class(&self) -> u16 {
        self.class
    }
    /// Returns the [`Chirality`] of the atom
    #[must_use]
    pub fn chirality(&self) -> Option<Chirality> {
        self.chirality
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
    #[must_use]
    pub fn with_isotope(mut self, iso: u16) -> Self {
        self.bracket_atom.isotope_mass_number = Some(iso);
        self
    }
    /// Adds the Atom Symbol
    #[must_use]
    pub fn with_symbol(mut self, symbol: AtomSymbol) -> Self {
        self.bracket_atom.symbol = symbol;
        self
    }
    /// Adds the aromatic value
    #[must_use]
    pub fn with_aromatic(mut self, aromatic: bool) -> Self {
        self.bracket_atom.aromatic = aromatic;
        self
    }
    /// Adds a specified [`HydrogenCount`]
    #[must_use]
    pub fn with_hydrogens(mut self, h_count: HydrogenCount) -> Self {
        self.bracket_atom.hydrogens = h_count;
        self
    }
    /// Adds a specified [`Charge`]
    #[must_use]
    pub fn with_charge(mut self, charge: Charge) -> Self {
        self.bracket_atom.charge = charge;
        self
    }
    /// Adds a designated esoteric class
    #[must_use]
    pub fn with_class(mut self, class: u16) -> Self {
        self.bracket_atom.class = class;
        self
    }
    /// Adds a specified [`Chirality`]
    #[must_use]
    pub fn with_chirality(mut self, chirality: Chirality) -> Self {
        self.bracket_atom.chirality = Some(chirality);
        self
    }
    /// Returns the [`Element`] contained in builder
    #[must_use]
    pub fn element(&self) -> Option<Element> {
        self.bracket_atom.element()
    }
    /// Returns the [`AtomSymbol`]
    #[must_use]
    pub fn symbol(&self) -> AtomSymbol {
        self.bracket_atom.symbol
    }
    /// Consumes the builder and returns the completed [`BracketAtom`]
    #[must_use]
    pub fn build(self) -> BracketAtom {
        BracketAtom {
            symbol: self.bracket_atom.symbol,
            aromatic: self.bracket_atom.aromatic,
            isotope_mass_number: self.bracket_atom.isotope_mass_number,
            hydrogens: self.bracket_atom.hydrogens,
            charge: self.bracket_atom.charge,
            class: self.bracket_atom.class,
            chirality: self.bracket_atom.chirality,
        }
    }
}

impl fmt::Display for BracketAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("[")?;
        if let Some(isotope) = self.isotope_mass_number() {
            write!(f, "{isotope}")?;
        }
        match self.symbol() {
            AtomSymbol::Element(element) => {
                if self.aromatic() {
                    write!(f, "{}", element.to_string().to_ascii_lowercase())?;
                } else {
                    write!(f,"{element}")?;
                }
            },
            AtomSymbol::WildCard | AtomSymbol::Unspecified => f.write_str("*")?,
        }
        if let Some(chirality) = self.chirality() {
            write!(f, "{chirality}")?;
        }
        write!(f, "{}", self.hydrogens())?;
        write!(f, "{}", self.charge())?;
        if self.class() != 0 {
            write!(f, ":{}", self.class())?;
        }
        f.write_str("]")
    }
}

#[cfg(test)]
mod tests {
    use elements_rs::Element;

    use super::BracketAtom;
    use crate::{
        atom::{
            atom_symbol::AtomSymbol,
            bracketed::{charge::Charge, chirality::Chirality, hydrogen_count::HydrogenCount},
        },
        errors::SmilesError,
    };

    #[test]
    fn builder_defaults_are_correct() {
        let a = BracketAtom::builder().build();

        assert_eq!(a.symbol(), AtomSymbol::default());
        assert_eq!(a.element(), None);
        assert_eq!(a.isotope_mass_number(), None);
        assert!(!a.aromatic());
        assert_eq!(a.hydrogens(), HydrogenCount::Unspecified);
        assert_eq!(a.hydrogen_count(), None);
        assert_eq!(a.charge_value(), 0);
        assert_eq!(a.class(), 0);
        assert_eq!(a.chirality(), None);
    }

    #[test]
    fn builder_setters_roundtrip_into_built_atom() {
        let sym = AtomSymbol::Element(Element::C);
        let charge = Charge::try_new(-3).unwrap();
        let chiral = Chirality::AtAt;

        let a = BracketAtom::builder()
            .with_symbol(sym)
            .with_isotope(13)
            .with_aromatic(true)
            .with_hydrogens(HydrogenCount::Explicit(2))
            .with_charge(charge)
            .with_class(42)
            .with_chirality(chiral)
            .build();

        assert_eq!(a.symbol(), sym);
        assert_eq!(a.element(), Some(Element::C));
        assert_eq!(a.isotope_mass_number(), Some(13));
        assert!(a.aromatic());
        assert_eq!(a.hydrogens(), HydrogenCount::Explicit(2));
        assert_eq!(a.hydrogen_count(), Some(2));
        assert_eq!(a.charge(), charge);
        assert_eq!(a.charge_value(), -3);
        assert_eq!(a.class(), 42);
        assert_eq!(a.chirality(), Some(chiral));
    }

    #[test]
    fn builder_element_and_symbol_accessors_reflect_current_state() {
        let b = BracketAtom::builder().with_symbol(AtomSymbol::Element(Element::N));

        assert_eq!(b.symbol(), AtomSymbol::Element(Element::N));
        assert_eq!(b.element(), Some(Element::N));
    }

    #[test]
    fn element_is_none_for_wildcard_and_unspecified() {
        let a = BracketAtom::builder().with_symbol(AtomSymbol::WildCard).build();
        assert_eq!(a.element(), None);

        let b = BracketAtom::builder().with_symbol(AtomSymbol::Unspecified).build();
        assert_eq!(b.element(), None);
    }

    #[test]
    fn isotope_returns_most_abundant_when_mass_is_none() {
        let a = BracketAtom::builder().with_symbol(AtomSymbol::Element(Element::C)).build();

        let iso = a.isotope().unwrap();
        assert_eq!(iso, Element::C.most_abundant_isotope());
    }

    #[test]
    fn isotope_returns_specific_when_mass_is_set() {
        let a = BracketAtom::builder()
            .with_symbol(AtomSymbol::Element(Element::C))
            .with_isotope(13)
            .build();

        let iso = a.isotope().unwrap();
        assert_eq!(iso, elements_rs::Isotope::try_from((Element::C, 13u16)).unwrap());
    }

    #[test]
    fn isotope_errors_when_symbol_has_no_element() {
        let a = BracketAtom::builder().with_symbol(AtomSymbol::WildCard).with_isotope(13).build();
        assert_eq!(a.isotope(), Err(SmilesError::InvalidIsotope));

        let b = BracketAtom::builder().with_symbol(AtomSymbol::Unspecified).build();
        assert_eq!(b.isotope(), Err(SmilesError::InvalidIsotope));
    }

    #[test]
    fn hydrogen_count_matches_hydrogens_enum() {
        let a = BracketAtom::builder().with_hydrogens(HydrogenCount::Unspecified).build();
        assert_eq!(a.hydrogen_count(), None);

        let b = BracketAtom::builder().with_hydrogens(HydrogenCount::Explicit(0)).build();
        assert_eq!(b.hydrogen_count(), Some(0));

        let c = BracketAtom::builder().with_hydrogens(HydrogenCount::Explicit(4)).build();
        assert_eq!(c.hydrogen_count(), Some(4));
    }

    #[test]
    fn charge_roundtrips_and_defaults_to_zero() {
        let a = BracketAtom::builder().build();
        assert_eq!(a.charge_value(), 0);

        let b = BracketAtom::builder().with_charge(Charge::try_new(5).unwrap()).build();
        assert_eq!(b.charge_value(), 5);
    }

    #[test]
    fn chiral_is_none_by_default_and_some_when_set() {
        let a = BracketAtom::builder().build();
        assert_eq!(a.chirality(), None);

        let b = BracketAtom::builder().with_chirality(Chirality::At).build();
        assert_eq!(b.chirality(), Some(Chirality::At));
    }

    #[test]
    fn build_consumes_builder_and_preserves_all_fields() {
        let a = BracketAtom::builder()
            .with_symbol(AtomSymbol::Element(Element::O))
            .with_isotope(18)
            .with_aromatic(false)
            .with_hydrogens(HydrogenCount::Explicit(1))
            .with_charge(Charge::try_new(-1).unwrap())
            .with_class(7)
            .with_chirality(Chirality::AtAt)
            .build();

        assert_eq!(a.element(), Some(Element::O));
        assert_eq!(a.isotope_mass_number(), Some(18));
        assert_eq!(a.hydrogen_count(), Some(1));
        assert_eq!(a.charge_value(), -1);
        assert_eq!(a.class(), 7);
        assert_eq!(a.chirality(), Some(Chirality::AtAt));
    }

#[test]
fn test_bracketed_atom_fmt_all_arms() {
    let cases = [
        (
            BracketAtom::builder().build(),
            "[*]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::WildCard)
                .build(),
            "[*]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .build(),
            "[C]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_aromatic(true)
                .build(),
            "[c]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_isotope(13)
                .build(),
            "[13C]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_chirality(Chirality::At)
                .build(),
            "[C@]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_chirality(Chirality::AtAt)
                .build(),
            "[C@@]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_hydrogens(HydrogenCount::Explicit(2))
                .build(),
            "[CH2]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_charge(Charge::try_new(1).unwrap())
                .build(),
            "[C+]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_charge(Charge::try_new(3).unwrap())
                .build(),
            "[C+3]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_charge(Charge::try_new(-1).unwrap())
                .build(),
            "[C-]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_charge(Charge::try_new(-3).unwrap())
                .build(),
            "[C-3]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_class(7)
                .build(),
            "[C:7]",
        ),
        (
            BracketAtom::builder()
                .with_symbol(AtomSymbol::Element(Element::C))
                .with_isotope(13)
                .with_aromatic(true)
                .with_chirality(Chirality::At)
                .with_hydrogens(HydrogenCount::Explicit(2))
                .with_charge(Charge::try_new(-1).unwrap())
                .with_class(7)
                .build(),
            "[13c@H2-:7]",
        ),
    ];

    for (atom, expected) in cases {
        assert_eq!(expected, atom.to_string());
    }
}
}
