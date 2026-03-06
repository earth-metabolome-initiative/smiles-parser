//! Second pass that parses the [`TokenWithSpan`]

use std::collections::HashMap;

use crate::{
    atom::{Atom, atom_node::AtomNode},
    bond::{Bond, ring_num::RingNum},
    errors::{SmilesError, SmilesErrorWithSpan},
    smiles::Smiles,
    token::{Token, TokenWithSpan},
};

/// Contains the slice of tokens being iterated on and current position in that
/// slice
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
    ///   is no current atom to branch from.
    ///
    /// - [`SmilesError::UnexpectedRightParentheses`]: Encountered `)` without a
    ///   matching `(`.
    ///
    /// - [`SmilesError::UnclosedBranch`]: A boundary was encountered while a
    ///   branch is still open.
    ///
    /// - [`SmilesError::InvalidRingNumber`]: A ring closure token was
    ///   encountered without a current atom or a ring closure is left unmatched
    ///   by end-of-input.
    ///
    /// - [`SmilesError::UnclosedRing`]: A boundary was encountered while there
    ///   are still open ring closures.
    ///
    /// - [`SmilesError::IncompleteBond`]: A bond token was parsed but no atom
    ///   proceeded to complete the bond.
    ///
    /// - Any other error will be emitted as a  [`SmilesErrorWithSpan`]
    pub fn parse(mut self) -> Result<Smiles, SmilesErrorWithSpan> {
        let mut smiles = Smiles::new();

        let mut next_id: usize = 0;
        let mut last_atom: Option<usize> = None;
        let mut pending_bond: Option<Bond> = None;
        let mut branch_stack: Vec<usize> = Vec::new();
        let mut ring_open: HashMap<RingNum, (usize, Option<Bond>)> = HashMap::new();
        let mut last_span: (usize, usize) = (0, 0);

        while let Some(token_with_span) = self.current() {
            let start = token_with_span.start();
            let end = token_with_span.end();
            last_span = (start, end);
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
                    branch_stack.push(anchor);
                }
                Token::NonBond => {
                    validate_non_bond(
                        pending_bond,
                        branch_stack.is_empty(),
                        ring_open.is_empty(),
                        last_span,
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
        parse_end_check(pending_bond, branch_stack.is_empty(), ring_open.is_empty(), last_span)?;

        Ok(smiles)
    }
}

fn parse_end_check(
    pending_bond: Option<Bond>,
    branch_stack_empty: bool,
    ring_open_empty: bool,
    last_span: (usize, usize),
) -> Result<(), SmilesErrorWithSpan> {
    let (start, end) = last_span;
    let start = start.min(end.saturating_sub(1));
    let end = end.max(start.saturating_add(1));

    if let Some(bond) = pending_bond {
        return Err(SmilesErrorWithSpan::new(SmilesError::IncompleteBond(bond), start, end));
    }
    if !branch_stack_empty {
        return Err(SmilesErrorWithSpan::new(SmilesError::UnclosedBranch, start, end));
    }
    if !ring_open_empty {
        return Err(SmilesErrorWithSpan::new(SmilesError::UnclosedRing, start, end));
    }
    Ok(())
}

fn validate_non_bond(
    pending_bond: Option<Bond>,
    branch_stack_empty: bool,
    ring_open_empty: bool,
    last_span: (usize, usize),
) -> Result<(), SmilesErrorWithSpan> {
    let (start, end) = last_span;
    let start = start.min(end.saturating_sub(1));
    let end = end.max(start.saturating_add(1));

    if let Some(bond) = pending_bond {
        return Err(SmilesErrorWithSpan::new(SmilesError::IncompleteBond(bond), start, end));
    }
    if !branch_stack_empty {
        return Err(SmilesErrorWithSpan::new(SmilesError::UnclosedBranch, start, end));
    }
    if !ring_open_empty {
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
mod tests {
    use elements_rs::Element;

    use crate::{
        atom::{atom_symbol::AtomSymbol, unbracketed::UnbracketedAtom},
        bond::{Bond, ring_num::RingNum},
        parser::smiles_parser::SmilesParser, // adjust if your module path differs
        token::{Token, TokenWithSpan},
    };

    fn test_tokens() -> Vec<TokenWithSpan> {
        let c = UnbracketedAtom::new(AtomSymbol::Element(Element::C), false);
        let o = UnbracketedAtom::new(AtomSymbol::Element(Element::O), false);
        let n = UnbracketedAtom::new(AtomSymbol::Element(Element::N), false);
        let aromatic_c = UnbracketedAtom::new(AtomSymbol::Element(Element::C), true);

        vec![
            // COc(c1)cccc1C#N
            TokenWithSpan::new(Token::UnbracketedAtom(c), 0, 1), // C
            TokenWithSpan::new(Token::UnbracketedAtom(o), 1, 2), // O
            TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 2, 3), // c
            TokenWithSpan::new(Token::LeftParentheses, 3, 4),    // (
            TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 4, 5), // c
            TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1).unwrap()), 5, 6), // 1
            TokenWithSpan::new(Token::RightParentheses, 6, 7),   // )
            TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 7, 8), // c
            TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 8, 9), // c
            TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 9, 10), // c
            TokenWithSpan::new(Token::UnbracketedAtom(aromatic_c), 10, 11), // c
            TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1).unwrap()), 11, 12), // 1
            TokenWithSpan::new(Token::UnbracketedAtom(c), 12, 13), // C
            TokenWithSpan::new(Token::Bond(Bond::Triple), 13, 14), // #
            TokenWithSpan::new(Token::UnbracketedAtom(n), 14, 15), // N
        ]
    }

    #[test]
    fn smiles_parser_new_sets_position_zero() {
        let tokens = test_tokens();
        let p = SmilesParser::new(&tokens);
        assert_eq!(p.position(), 0);
    }

    #[test]
    fn smiles_parser_tokens_returns_same_slice() {
        let tokens = test_tokens();
        let p = SmilesParser::new(&tokens);
        assert_eq!(p.tokens(), tokens.as_slice());
    }

    #[test]
    fn smiles_parser_current_peek_methods_work() {
        let tokens = test_tokens();
        let p = SmilesParser::new(&tokens);

        assert_eq!(p.current(), Some(&tokens[0]));
        assert_eq!(p.peek_next(), Some(&tokens[1]));
        assert_eq!(p.peek_n(2), Some(&tokens[2]));
        assert_eq!(p.peek_last(), None);
    }

    #[test]
    fn smiles_parser_advance_moves_position() {
        let tokens = test_tokens();
        let mut p = SmilesParser::new(&tokens);

        p.advance();
        assert_eq!(p.position(), 1);
        assert_eq!(p.current(), Some(&tokens[1]));
    }

    #[test]
    fn smiles_parser_next_token_returns_and_advances() {
        let tokens = test_tokens();
        let mut p = SmilesParser::new(&tokens);

        let t0 = p.next_token();
        assert_eq!(t0, Some(&tokens[0]));
        assert_eq!(p.position(), 1);

        let t1 = p.next_token();
        assert_eq!(t1, Some(&tokens[1]));
        assert_eq!(p.position(), 2);
    }

    #[test]
    fn smiles_parser_done_is_false_then_true() {
        let tokens = test_tokens();
        let mut p = SmilesParser::new(&tokens);

        assert!(!p.done());
        while !p.done() {
            p.advance();
        }
        assert!(p.done());
        assert_eq!(p.current(), None);
        assert_eq!(p.peek_next(), None);
        assert_eq!(p.peek_n(1), None);
    }

    #[test]
    fn smiles_parser_parse_builds_graph() {
        let tokens = test_tokens();
        let smiles = SmilesParser::new(&tokens).parse().expect("parse should succeed");

        assert_eq!(smiles.nodes().len(), 10);
        assert!(!smiles.edges().is_empty());
    }
}
