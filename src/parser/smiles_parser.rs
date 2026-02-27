//! Second pass that parses the [`TokenWithSpan`]

use crate::{
    atom::{Atom, atom_node::AtomNode, bracketed::BracketAtom, unbracketed::UnbracketedAtom}, bond::bond_edge::{self, BondEdge}, errors::{SmilesError, SmilesErrorWithSpan}, smiles::Smiles, token::{self, Token, TokenWithSpan}
};

/// Contains the vec of tokens being iterated on and tracks the current position
/// in that vec
pub struct SmilesParser<'a> {
    tokens: &'a [TokenWithSpan],
    position: usize,
}

impl<'a> SmilesParser<'a> {
    /// Creates a new `SmilesParser` structure
    #[must_use]
    pub fn new(tokens: &'a [TokenWithSpan]) -> Self {
        SmilesParser { tokens, position: 0 }
    }
    /// Retrieves the `tokens` field of [`Vec<TokenWithSpan>`]
    #[must_use]
    pub fn tokens(&self) -> &[TokenWithSpan] {
        self.tokens
    }
    /// Retrieves the current position
    #[must_use]
    pub fn position(&self) -> usize {
        self.position
    }
    /// Returns the current token, returns None if no token
    #[must_use]
    pub fn current(&self) -> Option<&TokenWithSpan> {
        self.tokens.get(self.position)
    }
    /// Returns the token `n` positions ahead of the current one without
    /// advancing.
    #[must_use]
    pub fn peek_n(&self, n: usize) -> Option<&TokenWithSpan> {
        self.tokens.get(self.position + n)
    }
    /// Returns the next token after the current one without advancing.
    #[must_use]
    pub fn peek_next(&self) -> Option<&TokenWithSpan> {
        self.tokens.get(self.position + 1)
    }
    /// Returns the last token before the current one without changing position
    #[must_use]
    pub fn peek_last(&self) -> Option<&TokenWithSpan> {
        self.position.checked_sub(1).and_then(|i| self.tokens.get(i))
    }
    /// Consumes and returns the next token
    pub fn next(&mut self) -> Option<&TokenWithSpan> {
        let token = self.tokens.get(self.position);
        if token.is_some() {
            self.position += 1;
        }
        token
    }
    /// Checks if token parsing has reached the end
    #[must_use]
    pub fn done(&self) -> bool {
        self.position >= self.tokens.len()
    }
    /// Advances the position
    pub fn advance(&mut self) {
        if !self.done() {
            self.position += 1;
        }
    }
    /// Parses the tokens to construct the [`Smiles`] structure
    pub fn parse(mut self) -> Result<Smiles, SmilesErrorWithSpan> {
        let mut smiles = Smiles::new();
        let mut id_count: usize = 0;
        while let Some(token_with_span) = self.current() {
            match token_with_span.token() {
                Token::NonBond => try_non_bond(self.peek_last(), self.peek_next(), token_with_span)?,
                Token::BracketedAtom(atom) => {
                    let atom = Atom::from(atom);

                    id_count +=1;
                },
                Token::UnbracketedAtom(atom) => todo!(),
                Token::Bond(bond) => todo!(),
                Token::LeftParentheses => todo!(),
                Token::RightParentheses => todo!(),
                Token::RingClosure(ring_num) => todo!(),
            }
            self.advance();
        }
        Ok(smiles)
    }

}

fn try_non_bond(
    prev_token: Option<&TokenWithSpan>,
    next_token: Option<&TokenWithSpan>,
    dot_token: &TokenWithSpan,
) -> Result<(), SmilesErrorWithSpan> {
    let prev_is_atom = prev_token.is_some_and(|t| {
        matches!(
            t.token(),
            Token::BracketedAtom(_) | Token::UnbracketedAtom(_)
        )
    });

    let next_is_atom = next_token.is_some_and(|t| {
        matches!(
            t.token(),
            Token::BracketedAtom(_) | Token::UnbracketedAtom(_)
        )
    });

    if !prev_is_atom || !next_is_atom {
        return Err(SmilesErrorWithSpan::new(
            SmilesError::InvalidNonBondToken,
            dot_token.start(),
            dot_token.end(),
        ));
    }

    Ok(())
}
