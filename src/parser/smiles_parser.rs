//! Second pass that parses the [`TokenWithSpan`]

use std::collections::HashMap;

use crate::{
    atom::{Atom, atom_node::AtomNode},
    bond::{Bond, ring_num::RingNum},
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
    pub fn next_token(&mut self) -> Option<&TokenWithSpan> {
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
    ///
    /// # Errors
    /// - [`SmilesError::UnexpectedLeftParentheses`]: Encountered `(` when there
    ///   is no current atom to branch from (e.g. the input begins with `(`, or
    ///   appears after `.`).
    ///
    /// - [`SmilesError::UnexpectedRightParentheses`]: Encountered `)` without a
    ///   matching `(` (branch stack underflow).
    ///
    /// - [`SmilesError::UnclosedBranch`]: A component boundary `.` was
    ///   encountered while there are still open branches.
    ///
    /// - [`SmilesError::InvalidRingNumber`]: A ring closure token was
    ///   encountered without a current atom (e.g. input starts with a ring
    ///   number), or a ring closure is left unmatched by end-of-input.
    ///
    /// - [`SmilesError::UnclosedRing`]: A component boundary `.` was
    ///   encountered while there are still open ring closures.
    ///
    /// - [`SmilesError::IncompleteBond`]: A bond token (e.g. `-`, `=`, `#`,
    ///   `:`) was read but no subsequent atom or ring closure was available to
    ///   complete the bond (for example, a bond immediately followed by `.`).
    ///
    /// - Any error produced while adding edges to the graph (e.g. invalid node
    ///   ids), which is wrapped into [`SmilesErrorWithSpan`] using the span of
    ///   the token that attempted the insertion.
    pub fn parse(mut self) -> Result<Smiles, SmilesErrorWithSpan> {
        let mut smiles = Smiles::new();

        let mut next_id: usize = 0;
        let mut last_atom: Option<usize> = None;
        let mut pending_bond: Option<Bond> = None;
        let mut branch_stack: Vec<usize> = Vec::new();
        let mut ring_open: HashMap<RingNum, (usize, Option<Bond>)> = HashMap::new();

        while let Some(token_with_span) = self.current() {
            let start = token_with_span.start();
            let end = token_with_span.end();
            match token_with_span.token() {
                Token::UnbracketedAtom(atom) => {
                    let atom: Atom = Atom::Unbracketed(atom);
                    let id = next_id;
                    next_id += 1;

                    let node = AtomNode::new(atom, id, token_with_span.span(), None);
                    smiles.push_node(node);

                    if let Some(prev) = last_atom {
                        let bond = pending_bond.unwrap_or_else(|| default_bond(&smiles, prev, id));
                        smiles
                            .push_edge(prev, id, bond)
                            .map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;
                    }
                    last_atom = Some(id);
                    pending_bond = None;
                }
                Token::BracketedAtom(atom) => {
                    let atom: Atom = Atom::from(atom);
                    let id = next_id;
                    next_id += 1;

                    let node = AtomNode::new(atom, id, token_with_span.span(), None);
                    smiles.push_node(node);

                    if let Some(prev) = last_atom {
                        let bond = pending_bond.unwrap_or_else(|| default_bond(&smiles, prev, id));
                        smiles
                            .push_edge(prev, id, bond)
                            .map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;
                    }

                    last_atom = Some(id);
                    pending_bond = None;
                }
                Token::Bond(bond) => pending_bond = Some(bond),
                Token::LeftParentheses => {
                    let Some(anchor) = last_atom else {
                        return Err(SmilesErrorWithSpan::new(
                            SmilesError::UnexpectedLeftParentheses,
                            start,
                            end,
                        ));
                    };
                    last_atom = Some(anchor);
                }
                Token::NonBond => {
                    is_valid_non_bond(
                        pending_bond,
                        start,
                        end,
                        branch_stack.is_empty(),
                        ring_open.is_empty(),
                    )?;
                    last_atom = None;
                    pending_bond = None;
                }
                Token::RingClosure(ring_num) => {
                    let Some(current) = last_atom else {
                        return Err(SmilesErrorWithSpan::new(
                            SmilesError::InvalidRingNumber,
                            start,
                            end,
                        ));
                    };
                    if let Some((other, stored_bond)) = ring_open.remove(&ring_num) {
                        let bond = pending_bond
                            .or(stored_bond)
                            .unwrap_or_else(|| default_bond(&smiles, current, other));
                        smiles
                            .push_edge(current, other, bond)
                            .map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;
                        pending_bond = None;
                    } else {
                        ring_open.insert(ring_num, (current, pending_bond));
                        pending_bond = None;
                    }
                }
                Token::RightParentheses => {
                    let Some(anchor) = branch_stack.pop() else {
                        return Err(SmilesErrorWithSpan::new(
                            SmilesError::UnexpectedRightParentheses,
                            start,
                            end,
                        ));
                    };
                    last_atom = Some(anchor);
                }
            }
            self.advance();
        }

        if let Some((_ring_num, _)) = ring_open.into_iter().next() {
            return Err(SmilesErrorWithSpan::new(SmilesError::InvalidRingNumber, 0, 0));
        }
        Ok(smiles)
    }
}

fn is_valid_non_bond(
    pending_bond: Option<Bond>,
    start: usize,
    end: usize,
    branch_stack: bool,
    ring_open: bool,
) -> Result<(), SmilesErrorWithSpan> {
    if let Some(bond) = pending_bond {
        return Err(SmilesErrorWithSpan::new(SmilesError::IncompleteBond(bond), start - 1, end));
    }
    if !branch_stack {
        return Err(SmilesErrorWithSpan::new(SmilesError::UnclosedBranch, start, end));
    }
    if !ring_open {
        return Err(SmilesErrorWithSpan::new(SmilesError::UnclosedRing, start, end));
    }
    Ok(())
}

fn default_bond(smiles: &Smiles, id_a: usize, id_b: usize) -> Bond {
    let node_a = &smiles.nodes()[id_a];
    let node_b = &smiles.nodes()[id_b];
    if node_a.atom().aromatic() && node_b.atom().aromatic() { Bond::Aromatic } else { Bond::Single }
}

#[cfg(test)]
mod tests {}
