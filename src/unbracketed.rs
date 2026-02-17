//! Module for containing an organic element that occurs outside of brackets
//! `[]`: `B, C, N, O, P, S, F, Cl, Br, I, *`.
use elements_rs::Element;

use crate::atom_symbol::AtomSymbol;

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
    pub const fn new(symbol: AtomSymbol, aromatic: bool) -> Self {
        Self { symbol, aromatic }
    }
    /// Returns the [`AtomSymbol`] of the atom
    pub fn symbol(&self) -> AtomSymbol {
        self.symbol
    }
    /// Returns the [`Element`] or `None` if `WildCard`
    pub fn element(&self) -> Option<Element> {
        self.symbol.element()
    }
    /// Returns true of aromatic
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }
    /// Returns true if `AtomSymbol` is [`AtomSymbol::WildCard`]
    pub fn is_wildcard(&self) -> bool {
        self.symbol.is_wildcard()
    }
}
