//! Module for the symbols representing an element in a `SMILES` string
use core::{cmp::Ordering, fmt};

use elements_rs::Element;

#[derive(Copy, Default, Debug, PartialEq, Clone, Eq, Hash)]
/// Enum to allow for standard elements or the `WildCard` variant, represented
/// as `*`
pub enum AtomSymbol {
    /// The explicitly named [`Element`]
    Element(Element),
    /// `WildCard` variant, described [here](http://opensmiles.org/opensmiles.html#inatoms)
    #[default]
    WildCard,
}

impl AtomSymbol {
    /// Creates an atom symbol
    #[must_use]
    pub fn new(element_type: Option<Element>) -> Self {
        match element_type {
            Some(element) => AtomSymbol::Element(element),
            None => AtomSymbol::default(),
        }
    }
    /// creates an `AtomSymbol` set as `WildCard`
    #[must_use]
    pub fn new_wildcard() -> Self {
        Self::WildCard
    }
    /// Verifies whether the symbol present is a wildcard
    #[must_use]
    pub fn is_wildcard(&self) -> bool {
        matches!(self, AtomSymbol::WildCard)
    }
    /// Returns either the [`Element`] or `None` if wildcard
    #[must_use]
    pub fn element(&self) -> Option<Element> {
        match self {
            AtomSymbol::Element(e) => Some(*e),
            AtomSymbol::WildCard => None,
        }
    }
    /// Consumes the `AtomSymbol` and returns the [`Element`] or `None` if
    /// `WildCard`
    #[must_use]
    pub fn into_element(self) -> Option<Element> {
        match self {
            AtomSymbol::Element(e) => Some(e),
            AtomSymbol::WildCard => None,
        }
    }
}

impl fmt::Display for AtomSymbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Element(e) => f.write_str(e.symbol()),
            Self::WildCard => f.write_str("*"),
        }
    }
}

impl PartialOrd for AtomSymbol {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AtomSymbol {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::WildCard, Self::WildCard) => Ordering::Equal,
            (Self::WildCard, Self::Element(_)) => Ordering::Less,
            (Self::Element(_), Self::WildCard) => Ordering::Greater,
            (Self::Element(left), Self::Element(right)) => left.symbol().cmp(right.symbol()),
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;

    use elements_rs::Element;

    use crate::atom::atom_symbol::AtomSymbol;

    #[test]
    fn test_atom_symbols_all() {
        let hydrogen = Element::H;
        let hydro_symbol = AtomSymbol::new(Some(hydrogen));
        assert!(!hydro_symbol.is_wildcard());
        assert_eq!(hydro_symbol.element(), Some(hydrogen));
        let into_hydro = hydro_symbol.into_element();
        assert_eq!(into_hydro, Some(hydrogen));

        let default = AtomSymbol::default();
        assert_eq!(default, AtomSymbol::WildCard);

        let wild = AtomSymbol::new_wildcard();
        assert!(wild.is_wildcard());
    }

    #[test]
    fn test_atom_symbol_fmt_all_arms() {
        let cases = [
            (AtomSymbol::Element(Element::H), "H"),
            (AtomSymbol::WildCard, "*"),
            (AtomSymbol::new(None), "*"),
        ];

        for (symbol, expected) in cases {
            assert_eq!(expected, symbol.to_string());
        }
    }

    #[test]
    fn atom_symbols_have_a_stable_total_order() {
        assert!(AtomSymbol::WildCard < AtomSymbol::Element(Element::C));
        assert!(AtomSymbol::Element(Element::C) < AtomSymbol::Element(Element::O));
    }
}
