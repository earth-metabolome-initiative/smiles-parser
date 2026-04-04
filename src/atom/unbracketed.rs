//! Module for containing an organic element that occurs outside of brackets
//! `[]`: `B, C, N, O, P, S, F, Cl, Br, I, *`.
use alloc::string::ToString;
use core::fmt;

use elements_rs::Element;

use crate::atom::atom_symbol::AtomSymbol;

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Structure for aliphatic atoms, aromatic or non aromatic
pub struct UnbracketedAtom {
    /// Unbracketed elements as [`Element`]
    symbol: AtomSymbol,
    /// Whether the atom is aromatic
    aromatic: bool,
}

impl UnbracketedAtom {
    /// Creates a new `UnbracketedAtom`
    #[must_use]
    pub const fn new(symbol: AtomSymbol, aromatic: bool) -> Self {
        Self { symbol, aromatic }
    }
    /// Returns the [`AtomSymbol`] of the atom
    #[must_use]
    pub fn symbol(&self) -> AtomSymbol {
        self.symbol
    }
    /// Returns the [`Element`] or `None` if `WildCard`
    #[must_use]
    pub fn element(&self) -> Option<Element> {
        self.symbol.element()
    }
    /// Returns true of aromatic
    #[must_use]
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }
    /// Returns true if `AtomSymbol` is [`AtomSymbol::WildCard`]
    #[must_use]
    pub fn is_wildcard(&self) -> bool {
        self.symbol.is_wildcard()
    }
}

impl fmt::Display for UnbracketedAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.symbol(), self.aromatic()) {
            (AtomSymbol::Element(element), true) => {
                write!(f, "{}", element.to_string().to_ascii_lowercase())
            }
            _ => write!(f, "{}", self.symbol()),
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;

    use elements_rs::Element;

    use crate::atom::{atom_symbol::AtomSymbol, unbracketed::UnbracketedAtom};

    fn symbols() -> &'static [AtomSymbol] {
        &[AtomSymbol::WildCard, AtomSymbol::Element(elements_rs::Element::Ac)]
    }
    #[test]
    fn test_all_unbracketed_files_and_impls() {
        let symbols = symbols();
        for symbol in symbols {
            let atom = UnbracketedAtom::new(*symbol, false);
            assert_eq!(symbol, &atom.symbol());
            assert_eq!(symbol.element(), atom.element());
            assert!(!atom.aromatic());
            if symbol == &AtomSymbol::WildCard {
                assert!(atom.is_wildcard());
            }
        }
    }

    #[test]
    fn test_unbracketed_atom_fmt_all_arms() {
        let cases = [
            (UnbracketedAtom::new(AtomSymbol::Element(Element::C), false), "C"),
            (UnbracketedAtom::new(AtomSymbol::Element(Element::C), true), "c"),
            (UnbracketedAtom::new(AtomSymbol::Element(Element::Cl), false), "Cl"),
            (UnbracketedAtom::new(AtomSymbol::Element(Element::Cl), true), "cl"),
            (UnbracketedAtom::new(AtomSymbol::WildCard, false), "*"),
            (UnbracketedAtom::new(AtomSymbol::WildCard, true), "*"),
        ];

        for (atom, expected) in cases {
            assert_eq!(expected, atom.to_string());
        }
    }
}
