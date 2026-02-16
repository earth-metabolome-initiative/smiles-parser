//! Module for the symbols representing an element in a `SMILES` string
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
