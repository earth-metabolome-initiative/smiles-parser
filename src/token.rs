//! Represents tokens used in parsing SMILES strings.

use std::ops::Range;

use crate::{
    atom::{bracketed::BracketAtom, unbracketed::UnbracketedAtom},
    bond::{Bond, ring_num::RingNum},
};

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// Represented with a `.`
    NonBond,
    /// Elements that occur inside of `[]`, structured as [`BracketAtom`]
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
        Self { token, span: start..end }
    }
    /// Returns the token
    #[must_use]
    pub fn token(&self) -> Token {
        self.token
    }
    /// Returns the span as [`Range`]
    #[must_use]
    pub fn span(&self) -> Range<usize> {
        self.span.start..self.span.end
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

#[cfg(test)]
mod tests {
    use elements_rs::Element;

    use super::{Token, TokenWithSpan};
    use crate::{
        atom::{atom_symbol::AtomSymbol, bracketed::BracketAtom, unbracketed::UnbracketedAtom},
        bond::{Bond, ring_num::RingNum},
    };

    #[test]
    fn token_variants_can_be_constructed_and_compared() {
        let bracket_atom =
            BracketAtom::builder().with_symbol(AtomSymbol::Element(Element::C)).build();
        let unbracketed_atom = UnbracketedAtom::new(AtomSymbol::Element(Element::O), false);
        let ring_num = RingNum::try_new(12).unwrap();

        let cases = [
            Token::NonBond,
            Token::BracketedAtom(bracket_atom),
            Token::UnbracketedAtom(unbracketed_atom),
            Token::Bond(Bond::Double),
            Token::LeftParentheses,
            Token::RightParentheses,
            Token::RingClosure(ring_num),
        ];

        assert_eq!(cases[0], Token::NonBond);
        assert_eq!(cases[1], Token::BracketedAtom(bracket_atom));
        assert_eq!(cases[2], Token::UnbracketedAtom(unbracketed_atom));
        assert_eq!(cases[3], Token::Bond(Bond::Double));
        assert_eq!(cases[4], Token::LeftParentheses);
        assert_eq!(cases[5], Token::RightParentheses);
        assert_eq!(cases[6], Token::RingClosure(ring_num));
    }

    #[test]
    fn token_with_span_new_and_accessors_work() {
        let token = Token::Bond(Bond::Triple);
        let token_with_span = TokenWithSpan::new(token, 3, 7);

        assert_eq!(token_with_span.token(), Token::Bond(Bond::Triple));
        assert_eq!(token_with_span.span(), 3..7);
        assert_eq!(token_with_span.start(), 3);
        assert_eq!(token_with_span.end(), 7);
    }

    #[test]
    fn token_with_span_preserves_complex_token_variants() {
        let bracket_atom =
            BracketAtom::builder().with_symbol(AtomSymbol::Element(Element::N)).build();
        let ring_num = RingNum::try_new(9).unwrap();

        let bracketed = TokenWithSpan::new(Token::BracketedAtom(bracket_atom), 0, 3);
        let ring = TokenWithSpan::new(Token::RingClosure(ring_num), 5, 6);
        let non_bond = TokenWithSpan::new(Token::NonBond, 10, 11);

        assert_eq!(bracketed.token(), Token::BracketedAtom(bracket_atom));
        assert_eq!(bracketed.span(), 0..3);

        assert_eq!(ring.token(), Token::RingClosure(ring_num));
        assert_eq!(ring.start(), 5);
        assert_eq!(ring.end(), 6);

        assert_eq!(non_bond.token(), Token::NonBond);
        assert_eq!(non_bond.span(), 10..11);
    }

    #[test]
    fn token_with_span_clone_and_eq_behave_as_expected() {
        let original = TokenWithSpan::new(Token::LeftParentheses, 2, 3);
        let cloned = original.clone();

        assert_eq!(original, cloned);
        assert_eq!(original.token(), Token::LeftParentheses);
        assert_eq!(cloned.span(), 2..3);
    }
}
