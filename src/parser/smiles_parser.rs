//! Second pass that parses the [`TokenWithSpan`]

use std::ops::Range;

use crate::{
    atom::{Atom, atom_node::AtomNode},
    bond::{self, Bond, ring_num::RingNum},
    errors::{SmilesError, SmilesErrorWithSpan},
    smiles::Smiles,
    token::{self, Token, TokenWithSpan},
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
        let mut first_node_id_in_sequence: Option<(usize, TokenWithSpan)> = None;
        let mut previous_node_id: Option<usize> = None;

        while let Some(token_with_span) = self.current().cloned() {
            match token_with_span.token() {
                Token::UnbracketedAtom(unbracketed_atom) => {
                    let atom = Atom::Unbracketed(unbracketed_atom);
                    let id = match previous_node_id {
                        Some(id) => id+1,
                        None => 0,
                    };
                    previous_node_id = Some(id);
                    let possible_ring: Option<RingNum> = if let Some(token) = self.peek_next() {
                        match token.token() {
                            Token::RingClosure(ring_num) => {self.advance(); Some(ring_num)},
                            _ => None,
                        }
                    } else {
                        None
                    };
                    let bond: Result<Option<Bond>, SmilesErrorWithSpan> = if let Some(token) = self.peek_next() {
                        match token.token(){
                            Token::NonBond => Ok(None),
                            Token::BracketedAtom(_) => Ok(Some(Bond::Single)),
                            Token::UnbracketedAtom(atom) => {
                                if unbracketed_atom.aromatic() && atom.aromatic() {
                                    Ok((Some(Bond::Aromatic)))
                                } else {
                                    Ok(Some(Bond::Single))
                                }
                            },
                            Token::Bond(bond) => {
                                self.advance();
                                Ok(Some(bond))},
                            Token::LeftParentheses => todo!(),
                            Token::RightParentheses => todo!(),
                            Token::RingClosure(ring_num) => Err(SmilesErrorWithSpan::new(SmilesError::InvalidRingNumber, token.start(), token.end())),
                        }
                    } else {
                        Ok(None)
                    };

                },
                Token::BracketedAtom(bracket_atom) => {
                    let atom = Atom::from(bracket_atom);
                    let current_node =
                        set_atom_node(previous_node_id, atom, None, token_with_span.span());
                    let current_id = current_node.id();
                    let current_val = (current_id, token_with_span);
                    first_node_id_in_sequence =
                        check_first_node(first_node_id_in_sequence, current_val);

                    previous_node_id = Some(current_id);
                    smiles.push_node(current_node);
                }
                Token::Bond(bond) => todo!(),
                Token::LeftParentheses => todo!(),
                Token::NonBond => todo!(),
                Token::RightParentheses => todo!(),
                Token::RingClosure(ring_num) => {}
            }
            self.advance();
        }

        Ok(smiles)
    }
}

fn set_atom_node(
    previous_node_id: Option<usize>,
    atom: Atom,
    ring_num: Option<RingNum>,
    span: Range<usize>,
) -> AtomNode {
    let id: usize;
    if let Some(prev) = previous_node_id {
        id = prev + 1;
    } else {
        id = 0;
    }
    let atom_node = AtomNode::new(atom, id, span, ring_num);
    atom_node
}

fn check_first_node(
    first_node_id_in_sequence: Option<(usize, TokenWithSpan)>,
    current_val: (usize, TokenWithSpan),
) -> Option<(usize, TokenWithSpan)> {
    if first_node_id_in_sequence.is_none() { Some(current_val) } else { first_node_id_in_sequence }
}

fn try_peek_bond(
    tokens: &[TokenWithSpan],
    next_token_location: usize,
) -> Result<(Bond, usize, usize), SmilesError> {
    let current = &tokens[next_token_location - 1];
    let next = &tokens[next_token_location];
    match next.token() {
        Token::NonBond => Ok((Bond::Single, next.start(), next.end())),
        Token::BracketedAtom(_) => Ok((Bond::Single, current.start(), next.end())),
        Token::UnbracketedAtom(_) => Ok((Bond::Single, current.start(), next.end())),
        Token::Bond(bond) => Ok((bond, next.start(), next.end())),
        Token::LeftParentheses => Ok((Bond::Single, current.start(), next.end())),
        Token::RightParentheses => Ok((Bond::Single, current.start(), next.end())),
        Token::RingClosure(ring_num) => {
            let second = &tokens[next_token_location + 1];
            match second.token() {
                Token::NonBond => todo!(),
                Token::BracketedAtom(bracket_atom) => todo!(),
                Token::UnbracketedAtom(unbracketed_atom) => todo!(),
                Token::Bond(bond) => todo!(),
                Token::LeftParentheses => todo!(),
                Token::RightParentheses => todo!(),
                Token::RingClosure(ring_num) => todo!(),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        errors::SmilesError, parser::smiles_parser::check_first_node, token::TokenWithSpan,
    };

    #[test]
    fn test_check_first_node() -> Result<(), SmilesError> {
        let node = 1;
        let token = TokenWithSpan::new(crate::token::Token::NonBond, 1, 2);
        let new_node = check_first_node(None, (node, token.clone()));
        if let Some((found_node, found_token)) = new_node {
            assert_eq!(found_node, node);
            assert_eq!(found_token, token);
        } else {
            return Err(SmilesError::NodeIdInvalid(node));
        }
        Ok(())
    }
}
