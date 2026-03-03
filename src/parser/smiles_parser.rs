//! Second pass that parses the [`TokenWithSpan`]

use std::collections::HashMap;

use crate::{
    atom::{Atom, atom_node::AtomNode},
    bond::Bond,
    errors::{SmilesError, SmilesErrorWithSpan},
    smiles::Smiles,
    token::{Token, TokenWithSpan},
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
        let mut first_node_id_in_sequence: Option<usize> = None;
        let mut previous_node_id: Option<usize> = None;
        // hold potential next bond and its span, needs to be evaluated as valid before pushing
        let mut pending_bond: Option<(Bond, usize, usize)> = None;

        while let Some(token_with_span) = self.current().cloned() {
            match token_with_span.token() {
                Token::NonBond => todo!(),
                Token::BracketedAtom(bracket_atom) => {
                    let atom = Atom::from(bracket_atom);
                    let current_node = set_atom_node(previous_node_id, atom);
                    check_first_node(&mut first_node_id_in_sequence, current_node.id());
                    previous_node_id = Some(current_node.id());
                    smiles.push_node(current_node);
                    

                },
                Token::UnbracketedAtom(unbracketed_atom) => todo!(),
                Token::Bond(bond) => todo!(),
                Token::LeftParentheses => todo!(),
                Token::RightParentheses => todo!(),
                Token::RingClosure(ring_num) => todo!(),
            }
            self.advance();
        }

        if let Some((_, start, end)) = pending_bond {
            return Err(SmilesErrorWithSpan::new(SmilesError::UnexpectedEndOfString, start, end));
        }

        Ok(smiles)
    }
    
    /// Tries to set the bond between the current atom and the previous atom
    /// 
    /// # Error
    /// - Will return [`SmilesError::IncompleteBond`] as a [`SmilesErrorWithSpan`]
    fn try_add_bond_edge(self, current_atom_id: usize, ) -> Result<(Bond, usize, usize), SmilesErrorWithSpan> {

    } 
}

fn set_atom_node(previous_node_id: Option<usize>, atom: Atom) -> AtomNode {
    let id: usize;
    if let Some(prev) = previous_node_id {
        id = prev + 1;
    } else {
        id = 0;
    }
    let atom_node = AtomNode::new(atom, id);
    atom_node
}

fn check_first_node(first_node_id_in_sequence: &mut Option<usize>, current_node_id: usize) {
    if first_node_id_in_sequence.is_none() {
        *first_node_id_in_sequence = Some(current_node_id);
    }
}



#[cfg(test)]
mod tests {
    use crate::parser::smiles_parser::check_first_node;

    #[test]
    fn test_check_first_node() {
        let mut node = None;
        check_first_node(&mut node, 3);
        assert_eq!(node, Some(3));
    }
}