//! Wrapper module for both bracketed and unbracketed atoms
pub mod atom_node;
pub mod atom_symbol;
pub mod bracketed;
pub mod unbracketed;

use core::fmt;

use elements_rs::{Element, Isotope};

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
            Self::Unbracketed(atom) => atom.aromatic(),
            Self::Bracketed(atom) => atom.aromatic(),
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
        self.charge().get()
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
            Self::Unbracketed(unbracketed_atom) => {
                let element = unbracketed_atom.element().ok_or(SmilesError::InvalidIsotope)?;
                Ok(element.most_abundant_isotope())
            }
            Self::Bracketed(bracket_atom) => bracket_atom.isotope(),
        }
    }
    /// Returns the element parsed for the atom or None if the atom is a
    /// wildcard (*)
    #[must_use]
    pub fn element(&self) -> Option<Element> {
        match self {
            Atom::Unbracketed(unbracketed_atom) => unbracketed_atom.element(),
            Atom::Bracketed(bracket_atom) => bracket_atom.element(),
        }
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unbracketed(atom) => write!(f, "{atom}"),
            Self::Bracketed(atom) => write!(f, "{atom}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;

    use elements_rs::Element;

    use super::*;

    fn ac_symbol() -> AtomSymbol {
        AtomSymbol::Element(Element::Ac)
    }

    #[test]
    fn from_unbracketed_sets_correct_variant_and_symbol_delegates() {
        let unbracketed = UnbracketedAtom::new(ac_symbol(), false);
        let atom: Atom = Atom::from(unbracketed);

        assert!(matches!(atom, Atom::Unbracketed(_)));
        assert_eq!(atom.symbol().element(), Some(Element::Ac));
        assert!(!atom.aromatic());
    }

    #[test]
    fn from_bracketed_sets_correct_variant_and_symbol_delegates() {
        let bracketed = BracketAtom::builder().with_symbol(ac_symbol()).build();

        let atom: Atom = Atom::from(bracketed);

        assert!(matches!(atom, Atom::Bracketed(_)));
        assert_eq!(atom.symbol().element(), Some(Element::Ac));
    }

    #[test]
    fn unbracketed_bracket_only_fields_are_defaults() {
        let atom: Atom = UnbracketedAtom::new(ac_symbol(), true).into();

        assert!(atom.aromatic());
        assert_eq!(atom.symbol().element(), Some(Element::Ac));

        assert_eq!(atom.chirality(), None);
        assert_eq!(atom.class(), 0);

        assert_eq!(atom.charge(), Charge::default());
        assert_eq!(atom.charge_value(), Charge::default().get());

        assert_eq!(atom.hydrogen_count(), None);
        assert_eq!(atom.hydrogens(), HydrogenCount::Unspecified);
    }

    #[test]
    fn bracketed_defaults_match_bracket_atom_defaults() {
        let atom: Atom = BracketAtom::builder().with_symbol(ac_symbol()).build().into();

        assert_eq!(atom.symbol().element(), Some(Element::Ac));
        assert_eq!(atom.chirality(), None);
        assert_eq!(atom.charge(), Charge::default());
        assert_eq!(atom.charge_value(), 0);
        assert_eq!(atom.hydrogens(), HydrogenCount::Unspecified);
        assert_eq!(atom.hydrogen_count(), None);
    }

    #[test]
    fn bracketed_with_charge_is_reflected_in_wrapper() {
        let plus_one = Charge::try_new(1_i8).expect("charge +1 should be valid");

        let atom: Atom =
            BracketAtom::builder().with_symbol(ac_symbol()).with_charge(plus_one).build().into();

        assert_eq!(atom.charge_value(), 1);
        assert_eq!(atom.charge(), plus_one);
    }

    #[test]
    fn bracketed_with_hydrogens_is_reflected_in_wrapper() {
        let h1 = HydrogenCount::Explicit(1);

        let atom: Atom =
            BracketAtom::builder().with_symbol(ac_symbol()).with_hydrogens(h1).build().into();

        assert_eq!(atom.hydrogens(), h1);
        assert_eq!(atom.hydrogen_count(), Some(1));
    }

    #[test]
    fn bracketed_with_class_is_reflected_in_wrapper() {
        let atom: Atom =
            BracketAtom::builder().with_symbol(ac_symbol()).with_class(12).build().into();

        assert_eq!(atom.class(), 12);
    }

    #[test]
    fn bracketed_with_chirality_is_reflected_in_wrapper() {
        let chirality = Chirality::At;

        let atom: Atom = BracketAtom::builder()
            .with_symbol(ac_symbol())
            .with_chirality(chirality)
            .build()
            .into();

        assert_eq!(atom.chirality(), Some(chirality));
    }

    #[test]
    fn isotope_unbracketed_ok_for_element() {
        let atom: Atom = UnbracketedAtom::new(ac_symbol(), false).into();
        assert!(atom.isotope().is_ok());
    }

    #[test]
    fn isotope_bracketed_delegates_to_bracket_atom() {
        let atom: Atom = BracketAtom::builder().with_symbol(ac_symbol()).build().into();

        assert!(atom.isotope().is_ok());
    }

    #[test]
    fn test_atom_fmt_all_arms() {
        let cases = [
            (Atom::from(UnbracketedAtom::new(AtomSymbol::Element(Element::C), false)), "C"),
            (Atom::from(UnbracketedAtom::new(AtomSymbol::Element(Element::C), true)), "c"),
            (Atom::from(UnbracketedAtom::new(AtomSymbol::WildCard, false)), "*"),
            (
                Atom::from(
                    BracketAtom::builder().with_symbol(AtomSymbol::Element(Element::C)).build(),
                ),
                "[C]",
            ),
            (
                Atom::from(
                    BracketAtom::builder()
                        .with_symbol(AtomSymbol::Element(Element::C))
                        .with_isotope(13)
                        .with_hydrogens(HydrogenCount::Explicit(2))
                        .with_charge(Charge::try_new(-1).unwrap())
                        .build(),
                ),
                "[13CH2-]",
            ),
        ];

        for (atom, expected) in cases {
            assert_eq!(expected, atom.to_string());
        }
    }
}
