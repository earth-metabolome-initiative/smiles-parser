//! Represents tokens used in parsing SMILES strings.

use std::ops::Range;

use elements_rs::Element;

use crate::{
    atom_symbol::AtomSymbol, bond::Bond, bracketed::bracket_atom::BracketAtom, errors::SmilesError,
};

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// Represented with a `.`
    NonBond,
    /// Elements that occur inside of `[]`, structured as [`BracketedAtom`]
    BracketedAtom(BracketAtom),
    /// Aliphatic organic elements only as [`UnbracketedAtom`]
    UnbracketedAtom(UnbracketedAtom),
    /// The parsed bond. Single bonds are implicit between atoms if not
    /// explicated as a bond, e.g. `CC` would be `C-C`.
    Bond(Bond),
    /// Token for left parentheses, used for branching `(`
    LeftParentheses,
    /// Token for right parentheses, used for ending a branch `)`
    RightParentheses,
    /// Ring number markers occur outside of `[]` and may be of type `%` and
    /// `0-99`, and `%` may be omitted.
    RingClosure(RingNum),
}

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
/// A parsed token and its relative location in the string
pub struct TokenWithSpan {
    /// The parsed token
    token: Token,
    /// The starting to ending points of the token in the string
    span: Range<usize>,
}

impl TokenWithSpan {
    /// Generates a new token with a specified span
    pub fn new(token: Token, start: usize, end: usize) -> Self {
        Self { token, span: Range { start, end } }
    }
    /// Returns the token
    pub fn token(&self) -> Token {
        self.token
    }
    /// Returns the span as [`Range`]
    pub fn span(&self) -> &Range<usize> {
        &self.span
    }
    /// Returns the start of the span as [`usize`]
    pub fn start(&self) -> usize {
        self.span.start
    }
    /// Returns the end of the span as [`usize`]
    pub fn end(&self) -> usize {
        self.span.end
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a ring marker and implements tighter bounds for the minimal and
/// maximal value a ring marker can be
pub struct RingNum(u8);
impl RingNum {
    /// Attempts to generate a [`RingNum`] form a [`u8`],
    ///
    /// # Errors
    /// - Returns a [`SmilesError::RingNumberOverflow`] if the value is above
    ///   `99`
    pub fn try_new(num: u8) -> Result<Self, SmilesError> {
        (0..=99).contains(&num).then_some(Self(num)).ok_or(SmilesError::RingNumberOverflow(num))
    }

    /// Returns the value for the [`RingNum`]
    pub fn get(&self) -> u8 {
        self.0
    }
}

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
