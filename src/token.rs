//! Represents tokens used in parsing SMILES strings.

use std::ops::Range;

use crate::{
    bond::Bond, bracketed::bracket_atom::BracketAtom, ring_num::RingNum,
    unbracketed::UnbracketedAtom,
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
    #[must_use]
    pub fn new(token: Token, start: usize, end: usize) -> Self {
        Self { token, span: Range { start, end } }
    }
    /// Returns the token
    #[must_use]
    pub fn token(&self) -> Token {
        self.token
    }
    /// Returns the span as [`Range`]
    #[must_use]
    pub fn span(&self) -> &Range<usize> {
        &self.span
    }
    /// Returns the start of the span as [`usize`]
    #[must_use]
    pub fn start(&self) -> usize {
        self.span.start
    }
    /// Returns the end of the span as [`usize`]
    #[must_use]
    pub fn end(&self) -> usize {
        self.span.end
    }
}
