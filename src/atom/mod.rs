//! Atom storage and helpers.
pub mod atom_symbol;
pub mod bracketed;

use alloc::{borrow::Cow, string::String};
use core::fmt;

use elements_rs::{Element, Isotope};
use geometric_traits::traits::TypedNode;

use crate::{
    atom::{
        atom_symbol::AtomSymbol,
        bracketed::{charge::Charge, chirality::Chirality},
    },
    errors::SmilesError,
};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
/// Distinguishes between the two SMILES atom syntaxes this crate stores.
pub enum AtomSyntax {
    /// Organic-subset atom written without brackets, such as `C` or `c`.
    OrganicSubset,
    /// Bracket atom written explicitly inside `[]`.
    Bracket,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
/// A parsed atom together with the syntax form it originated from.
pub struct Atom {
    symbol: AtomSymbol,
    isotope_mass_number: Option<u16>,
    aromatic: bool,
    /// Explicit hydrogen count written in the source.
    ///
    /// For bracket atoms, omitted `H` means `0`.
    /// For organic-subset atoms, this is always `0`.
    hydrogens: u8,
    charge: Charge,
    class: u16,
    chirality: Option<Chirality>,
    syntax: AtomSyntax,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
/// MCES-facing atom label derived from a parsed [`Atom`].
///
/// This intentionally captures chemistry-facing identity for the default
/// labeled-MCES bridge while excluding parser spelling details such as bracket
/// syntax, explicit hydrogen notation, atom class, and chirality tokens.
pub struct McesAtomType {
    symbol: AtomSymbol,
    aromatic: bool,
    isotope_mass_number: Option<u16>,
    formal_charge: i8,
}

impl McesAtomType {
    /// Returns the atom symbol used by the MCES label.
    #[inline]
    #[must_use]
    pub fn symbol(self) -> AtomSymbol {
        self.symbol
    }

    /// Returns whether the MCES label treats the atom as aromatic.
    #[inline]
    #[must_use]
    pub fn aromatic(self) -> bool {
        self.aromatic
    }

    /// Returns the isotope mass number, if one is part of the MCES label.
    #[inline]
    #[must_use]
    pub fn isotope_mass_number(self) -> Option<u16> {
        self.isotope_mass_number
    }

    /// Returns the formal charge used by the MCES label.
    #[inline]
    #[must_use]
    pub fn formal_charge(self) -> i8 {
        self.formal_charge
    }
}

impl Atom {
    /// Returns a builder for bracket atoms.
    #[inline]
    #[must_use]
    pub fn builder() -> AtomBuilder {
        AtomBuilder {
            atom: Self {
                symbol: AtomSymbol::default(),
                isotope_mass_number: None,
                aromatic: false,
                hydrogens: 0,
                charge: Charge::default(),
                class: 0,
                chirality: None,
                syntax: AtomSyntax::Bracket,
            },
        }
    }

    /// Creates a new organic-subset atom written without brackets.
    #[inline]
    #[must_use]
    pub fn new_organic_subset(symbol: AtomSymbol, aromatic: bool) -> Self {
        Self {
            symbol,
            isotope_mass_number: None,
            aromatic,
            hydrogens: 0,
            charge: Charge::default(),
            class: 0,
            chirality: None,
            syntax: AtomSyntax::OrganicSubset,
        }
    }

    /// Creates a new bracket atom with all parsed fields set explicitly.
    #[inline]
    #[must_use]
    pub(crate) fn new_bracket(
        symbol: AtomSymbol,
        isotope_mass_number: Option<u16>,
        aromatic: bool,
        hydrogens: u8,
        charge: Charge,
        class: u16,
        chirality: Option<Chirality>,
    ) -> Self {
        Self {
            symbol,
            isotope_mass_number,
            aromatic,
            hydrogens,
            charge,
            class,
            chirality,
            syntax: AtomSyntax::Bracket,
        }
    }

    /// Returns the syntax category used to parse this atom.
    #[inline]
    #[must_use]
    pub fn syntax(&self) -> AtomSyntax {
        self.syntax
    }

    /// Returns whether this atom was parsed from bracket syntax.
    #[inline]
    #[must_use]
    pub fn is_bracket_atom(&self) -> bool {
        self.syntax == AtomSyntax::Bracket
    }

    /// Returns whether this atom was parsed as an organic-subset atom.
    #[inline]
    #[must_use]
    pub fn is_organic_subset_atom(&self) -> bool {
        self.syntax == AtomSyntax::OrganicSubset
    }

    /// Returns the parsed atom symbol.
    #[inline]
    #[must_use]
    pub fn symbol(&self) -> AtomSymbol {
        self.symbol
    }

    /// Returns the parsed element, or `None` for `*`.
    #[inline]
    #[must_use]
    pub fn element(&self) -> Option<Element> {
        match self.symbol {
            AtomSymbol::WildCard => None,
            AtomSymbol::Element(element) => Some(element),
        }
    }

    /// Returns the parsed isotope mass number, if present.
    #[inline]
    #[must_use]
    pub fn isotope_mass_number(&self) -> Option<u16> {
        self.isotope_mass_number
    }

    /// Returns the resolved isotope for the atom.
    ///
    /// # Errors
    /// - Returns [`SmilesError::InvalidIsotope`] if no element is available.
    /// - Propagates `elements-rs` isotope lookup errors.
    pub fn isotope(&self) -> Result<Isotope, SmilesError> {
        let element = self.element().ok_or(SmilesError::InvalidIsotope)?;
        let isotope = match self.isotope_mass_number {
            None => element.most_abundant_isotope(),
            Some(mass) => Isotope::try_from((element, mass))?,
        };
        Ok(isotope)
    }

    /// Returns whether the atom is aromatic.
    #[inline]
    #[must_use]
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }

    /// Returns the explicit hydrogen count written on the atom.
    #[inline]
    #[must_use]
    pub fn hydrogen_count(&self) -> u8 {
        self.hydrogens
    }

    /// Returns the formal charge.
    #[inline]
    #[must_use]
    pub fn charge(&self) -> Charge {
        self.charge
    }

    /// Returns the formal charge as an `i8`.
    #[inline]
    #[must_use]
    pub fn charge_value(&self) -> i8 {
        self.charge.get()
    }

    /// Returns the atom class.
    #[inline]
    #[must_use]
    pub fn class(&self) -> u16 {
        self.class
    }

    /// Returns the chirality tag, if present.
    #[inline]
    #[must_use]
    pub fn chirality(&self) -> Option<Chirality> {
        self.chirality
    }

    #[inline]
    #[must_use]
    pub(crate) fn rendered_string(&self) -> String {
        let mut rendered = String::with_capacity(self.rendered_len_hint());
        self.write_smiles(&mut rendered)
            .unwrap_or_else(|_| unreachable!("writing to String cannot fail"));
        rendered
    }

    #[inline]
    #[must_use]
    pub(crate) fn rendered_cow(&self) -> Cow<'static, str> {
        if let Some(rendered) = self.rendered_static() {
            Cow::Borrowed(rendered)
        } else {
            Cow::Owned(self.rendered_string())
        }
    }

    #[inline]
    #[must_use]
    pub(crate) fn rendered_len_hint(&self) -> usize {
        let mut len = rendered_symbol_len(self.symbol, self.aromatic);
        if self.syntax == AtomSyntax::Bracket {
            len += 2;
            if let Some(isotope) = self.isotope_mass_number {
                len += decimal_len_u16(isotope);
            }
            if let Some(chirality) = self.chirality {
                len += chirality.display_len();
            }
            if self.hydrogens != 0 {
                len += 1;
                if self.hydrogens != 1 {
                    len += decimal_len_u8(self.hydrogens);
                }
            }
            len += self.charge.display_len();
            if self.class != 0 {
                len += 1 + decimal_len_u16(self.class);
            }
        }
        len
    }

    pub(crate) fn write_smiles<W: fmt::Write>(&self, target: &mut W) -> fmt::Result {
        match self.syntax {
            AtomSyntax::OrganicSubset => write_symbol(target, self.symbol, self.aromatic),
            AtomSyntax::Bracket => {
                target.write_str("[")?;
                if let Some(isotope) = self.isotope_mass_number {
                    write!(target, "{isotope}")?;
                }
                write_symbol(target, self.symbol, self.aromatic)?;
                if let Some(chirality) = self.chirality {
                    write!(target, "{chirality}")?;
                }
                match self.hydrogens {
                    0 => {}
                    1 => target.write_str("H")?,
                    count => write!(target, "H{count}")?,
                }
                write!(target, "{}", self.charge)?;
                if self.class != 0 {
                    write!(target, ":{}", self.class)?;
                }
                target.write_str("]")
            }
        }
    }

    #[inline]
    fn rendered_static(&self) -> Option<&'static str> {
        match self.syntax {
            AtomSyntax::OrganicSubset => rendered_symbol_static(self.symbol, self.aromatic),
            AtomSyntax::Bracket => None,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
/// Builder for bracket atoms.
pub struct AtomBuilder {
    atom: Atom,
}

impl AtomBuilder {
    /// Adds an isotope value.
    #[inline]
    #[must_use]
    pub fn with_isotope(mut self, iso: u16) -> Self {
        self.atom.isotope_mass_number = Some(iso);
        self
    }

    /// Adds the atom symbol.
    #[inline]
    #[must_use]
    pub fn with_symbol(mut self, symbol: AtomSymbol) -> Self {
        self.atom.symbol = symbol;
        self
    }

    /// Adds aromaticity.
    #[inline]
    #[must_use]
    pub fn with_aromatic(mut self, aromatic: bool) -> Self {
        self.atom.aromatic = aromatic;
        self
    }

    /// Adds the explicit hydrogen count.
    #[inline]
    #[must_use]
    pub fn with_hydrogens(mut self, h_count: u8) -> Self {
        self.atom.hydrogens = h_count;
        self
    }

    /// Adds a formal charge.
    #[inline]
    #[must_use]
    pub fn with_charge(mut self, charge: Charge) -> Self {
        self.atom.charge = charge;
        self
    }

    /// Adds an atom class.
    #[inline]
    #[must_use]
    pub fn with_class(mut self, class: u16) -> Self {
        self.atom.class = class;
        self
    }

    /// Adds a chirality tag.
    #[inline]
    #[must_use]
    pub fn with_chirality(mut self, chirality: Chirality) -> Self {
        self.atom.chirality = Some(chirality);
        self
    }

    /// Returns the parsed element, or `None` for `*`.
    #[inline]
    #[must_use]
    pub fn element(&self) -> Option<Element> {
        self.atom.element()
    }

    /// Returns the current symbol.
    #[inline]
    #[must_use]
    pub fn symbol(&self) -> AtomSymbol {
        self.atom.symbol
    }

    /// Consumes the builder and returns the completed atom.
    #[inline]
    #[must_use]
    pub fn build(self) -> Atom {
        self.atom
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_smiles(f)
    }
}

impl TypedNode for Atom {
    type NodeType = McesAtomType;

    #[inline]
    fn node_type(&self) -> Self::NodeType {
        McesAtomType {
            symbol: self.symbol(),
            aromatic: self.aromatic(),
            isotope_mass_number: self.isotope_mass_number(),
            formal_charge: self.charge_value(),
        }
    }
}

#[inline]
fn write_symbol<W: fmt::Write>(target: &mut W, symbol: AtomSymbol, aromatic: bool) -> fmt::Result {
    match rendered_symbol_static(symbol, aromatic) {
        Some(rendered) => target.write_str(rendered),
        None => write!(target, "{symbol}"),
    }
}

#[inline]
fn rendered_symbol_len(symbol: AtomSymbol, aromatic: bool) -> usize {
    match rendered_symbol_static(symbol, aromatic) {
        Some(rendered) => rendered.len(),
        None => {
            match symbol {
                AtomSymbol::WildCard => 1,
                AtomSymbol::Element(element) => element.symbol_len(),
            }
        }
    }
}

#[inline]
fn rendered_symbol_static(symbol: AtomSymbol, aromatic: bool) -> Option<&'static str> {
    match (symbol, aromatic) {
        (AtomSymbol::WildCard, _) => Some("*"),
        (AtomSymbol::Element(element), false) => Some(element.symbol()),
        (AtomSymbol::Element(element), true) => element.aromatic_smiles_symbol(),
    }
}

const fn decimal_len_u8(value: u8) -> usize {
    if value >= 100 {
        3
    } else if value >= 10 {
        2
    } else {
        1
    }
}

const fn decimal_len_u16(value: u16) -> usize {
    if value >= 10_000 {
        5
    } else if value >= 1_000 {
        4
    } else if value >= 100 {
        3
    } else if value >= 10 {
        2
    } else {
        1
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;

    use elements_rs::Element;
    use geometric_traits::traits::TypedNode;

    use super::*;

    fn ac_symbol() -> AtomSymbol {
        AtomSymbol::Element(Element::Ac)
    }

    #[test]
    fn organic_subset_constructor_sets_defaults() {
        let atom = Atom::new_organic_subset(ac_symbol(), true);

        assert_eq!(atom.syntax(), AtomSyntax::OrganicSubset);
        assert!(atom.is_organic_subset_atom());
        assert!(!atom.is_bracket_atom());
        assert_eq!(atom.symbol().element(), Some(Element::Ac));
        assert!(atom.aromatic());
        assert_eq!(atom.chirality(), None);
        assert_eq!(atom.class(), 0);
        assert_eq!(atom.charge(), Charge::default());
        assert_eq!(atom.charge_value(), 0);
        assert_eq!(atom.hydrogen_count(), 0);
        assert_eq!(atom.isotope_mass_number(), None);
    }

    #[test]
    fn builder_defaults_are_bracket_atom_defaults() {
        let atom = Atom::builder().build();

        assert_eq!(atom.syntax(), AtomSyntax::Bracket);
        assert!(atom.is_bracket_atom());
        assert!(!atom.is_organic_subset_atom());
        assert_eq!(atom.symbol(), AtomSymbol::WildCard);
        assert!(!atom.aromatic());
        assert_eq!(atom.hydrogen_count(), 0);
        assert_eq!(atom.charge(), Charge::default());
        assert_eq!(atom.class(), 0);
        assert_eq!(atom.chirality(), None);
    }

    #[test]
    fn builder_setters_roundtrip_into_built_atom() {
        let atom = Atom::builder()
            .with_symbol(ac_symbol())
            .with_isotope(227)
            .with_aromatic(true)
            .with_hydrogens(4)
            .with_charge(Charge::try_new(3).unwrap())
            .with_class(12)
            .with_chirality(Chirality::At)
            .build();

        assert_eq!(atom.syntax(), AtomSyntax::Bracket);
        assert_eq!(atom.symbol().element(), Some(Element::Ac));
        assert_eq!(atom.isotope_mass_number(), Some(227));
        assert!(atom.aromatic());
        assert_eq!(atom.hydrogen_count(), 4);
        assert_eq!(atom.charge_value(), 3);
        assert_eq!(atom.class(), 12);
        assert_eq!(atom.chirality(), Some(Chirality::At));
    }

    #[test]
    fn builder_element_and_symbol_accessors_reflect_current_state() {
        let builder = Atom::builder().with_symbol(AtomSymbol::Element(Element::N));

        assert_eq!(builder.element(), Some(Element::N));
        assert_eq!(builder.symbol(), AtomSymbol::Element(Element::N));
    }

    #[test]
    fn element_is_none_for_wildcard_and_present_for_element() {
        let wildcard = Atom::builder().with_symbol(AtomSymbol::WildCard).build();
        assert_eq!(wildcard.element(), None);

        let element = Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).build();
        assert_eq!(element.element(), Some(Element::C));
    }

    #[test]
    fn isotope_returns_most_abundant_when_mass_is_none() {
        let atom = Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).build();
        assert_eq!(atom.isotope().unwrap(), Element::C.most_abundant_isotope());
    }

    #[test]
    fn isotope_returns_specific_when_mass_is_set() {
        let atom =
            Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).with_isotope(13).build();

        assert_eq!(atom.isotope().unwrap(), Isotope::try_from((Element::C, 13_u16)).unwrap());
    }

    #[test]
    fn isotope_errors_when_symbol_has_no_element() {
        let atom = Atom::builder().with_symbol(AtomSymbol::WildCard).build();
        assert_eq!(atom.isotope().unwrap_err(), SmilesError::InvalidIsotope);
    }

    #[test]
    fn typed_node_ignores_parser_spelling_fields() {
        let organic = Atom::new_organic_subset(AtomSymbol::Element(Element::C), false);
        let bracket = Atom::builder()
            .with_symbol(AtomSymbol::Element(Element::C))
            .with_hydrogens(3)
            .with_class(7)
            .with_chirality(Chirality::At)
            .build();

        assert_eq!(organic.node_type(), bracket.node_type());
    }

    #[test]
    fn typed_node_keeps_charge_isotope_and_aromaticity() {
        let neutral = Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).build();
        let charged = Atom::builder()
            .with_symbol(AtomSymbol::Element(Element::C))
            .with_charge(Charge::try_new(1).unwrap())
            .build();
        let isotopic =
            Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).with_isotope(13).build();
        let aromatic = Atom::new_organic_subset(AtomSymbol::Element(Element::C), true);
        let aliphatic = Atom::new_organic_subset(AtomSymbol::Element(Element::C), false);

        assert_ne!(neutral.node_type(), charged.node_type());
        assert_ne!(neutral.node_type(), isotopic.node_type());
        assert_ne!(aromatic.node_type(), aliphatic.node_type());
    }

    #[test]
    fn test_atom_fmt_all_arms() {
        let cases = [
            (Atom::new_organic_subset(AtomSymbol::Element(Element::C), false), "C"),
            (Atom::new_organic_subset(AtomSymbol::Element(Element::C), true), "c"),
            (Atom::new_organic_subset(AtomSymbol::WildCard, false), "*"),
            (Atom::builder().build(), "[*]"),
            (Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).build(), "[C]"),
            (
                Atom::builder()
                    .with_symbol(AtomSymbol::Element(Element::C))
                    .with_aromatic(true)
                    .build(),
                "[c]",
            ),
            (
                Atom::builder()
                    .with_symbol(AtomSymbol::Element(Element::C))
                    .with_hydrogens(1)
                    .build(),
                "[CH]",
            ),
            (
                Atom::builder()
                    .with_symbol(AtomSymbol::Element(Element::C))
                    .with_hydrogens(4)
                    .build(),
                "[CH4]",
            ),
            (
                Atom::builder()
                    .with_symbol(AtomSymbol::Element(Element::C))
                    .with_charge(Charge::try_new(2).unwrap())
                    .build(),
                "[C+2]",
            ),
            (
                Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).with_class(7).build(),
                "[C:7]",
            ),
            (
                Atom::builder()
                    .with_symbol(AtomSymbol::Element(Element::C))
                    .with_isotope(13)
                    .with_chirality(Chirality::AtAt)
                    .with_hydrogens(1)
                    .with_charge(Charge::try_new(-1).unwrap())
                    .with_class(2)
                    .build(),
                "[13C@@H-:2]",
            ),
        ];

        for (atom, expected) in cases {
            assert_eq!(atom.to_string(), expected);
        }
    }
}
