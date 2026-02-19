//! Module for the symbols representing an element in a `SMILES` string
use core::fmt;

use elements_rs::Element;

#[derive(Copy, Default, Debug, PartialEq, Clone, Eq, Hash)]
/// Enum to allow for standard elements or the `WildCard` variant, represented
/// as `*`
pub enum AtomSymbol {
    /// The explicitly named [`Element`]
    Element(Element),
    /// WildCard variant, described [here](http://opensmiles.org/opensmiles.html#inatoms)
    WildCard,
    /// unspecified Atom Symbol
    #[default]
    Unspecified,
}

impl AtomSymbol {
    /// Creates an atom symbol
    pub fn new(element_type: Option<Element>) -> Self {
        match element_type {
            Some(element) => AtomSymbol::Element(element),
            None => AtomSymbol::default(),
        }
    }
    /// creates an `AtomSymbol` set as `WildCard`
    pub fn new_wildcard() -> Self {
        Self::WildCard
    }
    /// Verifies whether the symbol present is a wildcard
    pub fn is_wildcard(&self) -> bool {
        matches!(self, AtomSymbol::WildCard)
    }
    /// Returns either the [`Element`] or `None` if wildcard
    pub fn element(&self) -> Option<Element> {
        match self {
            AtomSymbol::Element(e) => Some(*e),
            AtomSymbol::WildCard | AtomSymbol::Unspecified => None,
        }
    }
    /// Consumes the `AtomSymbol` and returns the [`Element`] or `None` if
    /// `WildCard`
    pub fn into_element(self) -> Option<Element> {
        match self {
            AtomSymbol::Element(e) => Some(e),
            AtomSymbol::WildCard | AtomSymbol::Unspecified => None,
        }
    }
}

impl fmt::Display for AtomSymbol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

#[cfg(test)]
mod tests {
    use elements_rs::Element;

    use crate::atom_symbol::AtomSymbol;

    #[test]
    fn test_atom_symbols_all() {
        let hydrogen = Element::H;
        let hydro_symbol = AtomSymbol::new(Some(hydrogen));
        assert!(!hydro_symbol.is_wildcard());
        assert_eq!(hydro_symbol.element(), Some(hydrogen));
        let into_hydro = hydro_symbol.into_element();
        assert_eq!(into_hydro, Some(hydrogen));

        let default = AtomSymbol::default();
        assert_eq!(default, AtomSymbol::Unspecified);

        let wild = AtomSymbol::new_wildcard();
        assert!(wild.is_wildcard());
    }
}
