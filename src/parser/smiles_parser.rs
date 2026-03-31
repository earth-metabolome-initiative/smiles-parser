//! Second pass that parses the [`TokenWithSpan`]

use std::{collections::HashMap, ops::Range};

use crate::{
    atom::{Atom, atom_node::AtomNode},
    bond::{Bond, ring_num::RingNum},
    errors::{SmilesError, SmilesErrorWithSpan},
    smiles::Smiles,
    token::{Token, TokenWithSpan},
};

/// Structure containing parser state.
pub struct ParserState {
    smiles: Smiles,
    next_id: usize,
    last_atom: Option<usize>,
    pending_bond: Option<Bond>,
    branch_stack: Vec<usize>,
    ring_open: HashMap<RingNum, (usize, Option<Bond>)>,
    last_span: (usize, usize),
    branch_status: (bool, bool),
}
impl ParserState {
    /// Creates a new initial state for the parser.
    #[must_use]
    pub fn new() -> Self {
        Self {
            smiles: Smiles::new(),
            next_id: 0,
            last_atom: None,
            pending_bond: None,
            branch_stack: Vec::new(),
            ring_open: HashMap::new(),
            last_span: (0, 0),
            branch_status: (false, false)
        }
    }
    /// Returns the last span stored.
    #[must_use]
    pub fn last_span(&self) -> (usize, usize) {
        self.last_span
    }
    /// Updates the last span field.
    pub fn update_last_span(&mut self, last_span: (usize, usize)) {
        self.last_span = last_span;
    }
    /// Returns the next ID field.
    #[must_use]
    pub fn next_id(&self) -> usize {
        self.next_id
    }
    /// Updates the next available ID by incrementing.
    pub fn increment_next_id(&mut self) {
        self.next_id += 1;
    }
    /// Returns the previous atom ID if there is one.
    #[must_use]
    pub fn last_atom(&self) -> Option<usize> {
        self.last_atom
    }
    /// Updates the last atom ID.
    pub fn update_last_atom(&mut self, id: Option<usize>) {
        self.last_atom = id;
    }
    /// Returns the pending bond if present
    #[must_use]
    pub fn pending_bond(&self) -> Option<Bond> {
        self.pending_bond
    }
    /// Updates the pending bond field.
    pub fn update_pending_bond(&mut self, bond: Option<Bond>) {
        self.pending_bond = bond;
    }
    /// Returns a borrowed slice of the current branch stack.
    #[must_use]
    pub fn branch_stack(&self) -> &[usize] {
        &self.branch_stack
    }
    /// Pops a value off the the branch stack and returns it.
    pub fn pop_branch_stack(&mut self) -> Option<usize> {
        self.branch_stack.pop()
    }
    /// Pushes a branch anchor to the branch stack.
    pub fn push_stack(&mut self, anchor: usize) {
        self.branch_stack.push(anchor);
    }
    /// Checks if the branch stack is empty.
    #[must_use]
    pub fn stack_empty(&self) -> bool {
        self.branch_stack.is_empty()
    }
    /// Removes and returns the specified ring open field entry if present.
    pub fn remove_ring_open(&mut self, ring_num: &RingNum) -> Option<(usize, Option<Bond>)> {
        self.ring_open.remove(ring_num)
    }
    /// Checks if the ring open field is currently empty.
    #[must_use]
    pub fn ring_open_empty(&self) -> bool {
        self.ring_open.is_empty()
    }
    /// Inserts the given ring into the ring open field
    pub fn insert_ring(&mut self, ring_num: RingNum, pending: (usize, Option<Bond>)) {
        self.ring_open.insert(ring_num, pending);
    }
    /// Returns a borrowed reference to the smiles field.
    #[must_use]
    pub fn smiles(&self) -> &Smiles {
        &self.smiles
    }
    /// Consumes the parser state and returns the parsed SMILES graph.
    #[must_use]
    pub fn into_smiles(self) -> Smiles {
        self.smiles
    }
    /// Returns whether there is an edge for the given pair of nodes.
    #[must_use]
    pub fn edge_for_node_pair_exists(&self, nodes: (usize, usize)) -> bool {
        self.smiles.edge_for_node_pair(nodes).is_some()
    }
    /// Pushes an [`AtomNode`] into the parsed [`Smiles`] graph.
    ///
    /// # Errors
    /// - Returns [`SmilesError::DuplicateNodeId`] if node ID already exists
    pub fn push_node(&mut self, node: AtomNode) -> Result<(), SmilesError> {
        self.smiles.push_node(node)
    }
    /// Pushes an edge into the parsed [`Smiles`] graph.
    ///
    /// # Errors
    /// - Returns a [`SmilesError::NodeIdInvalid`] if a node cannot be found in
    ///   the edge list
    pub fn push_edge(
        &mut self,
        node_a: usize,
        node_b: usize,
        bond: Bond,
        ring_num: Option<RingNum>,
    ) -> Result<(), SmilesError> {
        self.smiles.push_edge(node_a, node_b, bond, ring_num)
    }
    /// Adds an atom to the SMILES graph, either bracketed or unbracketed.
    ///
    /// # Errors
    /// - Returns [`SmilesError::DuplicateNodeId`] if node ID already exists
    pub fn add_atom(
        &mut self,
        atom: Atom,
        token_span: Range<usize>,
        start: usize,
        end: usize,
    ) -> Result<(), SmilesErrorWithSpan> {
        let id = self.next_id();
        self.increment_next_id();
        let node = AtomNode::new(atom, id, token_span);
        self.push_node(node).map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;
        if let Some(prev) = self.last_atom() {
            let bond = self.pending_bond().unwrap_or_else(|| default_bond(self.smiles(), prev, id));
            self.push_edge(prev, id, bond, None)
                .map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;
        }
        self.update_last_atom(Some(id));
        self.update_pending_bond(None);
        if self.branch_status.0 == true {
            self.branch_status.1 = true;
        }
        Ok(())
    }
    /// Validates that at the current point in parsing there are no hanging
    /// bonds, branches, or ring closures.
    ///
    /// # Errors
    /// - Returns [`SmilesError::IncompleteBond`] if a bond is left open at the
    ///   non-bond.
    /// - Returns [`SmilesError::UnclosedBranch`] if a branch is left open at
    ///   the non-bond.
    /// - Returns [`SmilesError::UnclosedRing`] if a ring is left open at the
    ///   non-bond.
    pub fn validate_all_closed(&mut self) -> Result<(), SmilesErrorWithSpan> {
        let (start, end) = self.last_span;
        let start = start.min(end.saturating_sub(1));
        let end = end.max(start.saturating_add(1));

        if let Some(bond) = self.pending_bond {
            return Err(SmilesErrorWithSpan::new(SmilesError::IncompleteBond(bond), start, end));
        }
        if !self.stack_empty() {
            return Err(SmilesErrorWithSpan::new(SmilesError::UnclosedBranch, start, end));
        }
        if !self.ring_open_empty() {
            return Err(SmilesErrorWithSpan::new(SmilesError::UnclosedRing, start, end));
        }
        self.update_last_atom(None);
        self.update_pending_bond(None);
        Ok(())
    }

    /// This method validates and adds the ring number to the relevant bond in
    /// the graph.
    ///
    /// # Errors
    /// - Returns [`SmilesError::InvalidRingNumber`] if a previous atom for the
    ///   bond is not found or a relevant edge between the vertices is not
    ///   found.
    /// - Returns [`SmilesError::NodeIdInvalid`] if a node cannot be found in
    ///   the edge list
    pub fn validate_and_add_ring_num(
        &mut self,
        start: usize,
        end: usize,
        ring_num: RingNum,
    ) -> Result<(), SmilesErrorWithSpan> {
        let Some(current) = self.last_atom() else {
            return Err(SmilesErrorWithSpan::new(SmilesError::InvalidRingNumber, start, end));
        };
        if let Some((other, stored_bond)) = self.remove_ring_open(&ring_num) {
            if current == other {
                return Err(SmilesErrorWithSpan::new(SmilesError::InvalidRingNumber, start, end));
            }
            if self.edge_for_node_pair_exists((current, other)) {
                return Err(SmilesErrorWithSpan::new(SmilesError::InvalidRingNumber, start, end));
            }
            let bond = self
                .pending_bond()
                .or(stored_bond)
                .unwrap_or_else(|| default_bond(self.smiles(), current, other));

            self.push_edge(current, other, bond, Some(ring_num))
                .map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;

            self.update_pending_bond(None);
        } else {
            self.insert_ring(ring_num, (current, self.pending_bond()));
            self.update_pending_bond(None);
        }

        Ok(())
    }
    /// Evaluates that a left token parentheses has a valid anchor atom.
    ///
    /// # Errors
    /// - Returns [`SmilesError::UnexpectedLeftParentheses`] if a valid anchor
    ///   is not found to associate with the left parentheses.
    pub fn validate_branch_open(
        &mut self,
        start: usize,
        end: usize,
        next_token: Option<&TokenWithSpan>,
    ) -> Result<(), SmilesErrorWithSpan> {
        self.branch_status.0 = true;
        if let Some(token) = next_token {
            if token.token() == Token::LeftParentheses {
                return Err(SmilesErrorWithSpan::new(
                    SmilesError::UnexpectedLeftParentheses,
                    start,
                    end,
                ));
            } else if token.token() == Token::RightParentheses {
                return Err(SmilesErrorWithSpan::new(SmilesError::EmptyBranch, start, end));
            }
        }
        let Some(anchor) = self.last_atom() else {
            return Err(SmilesErrorWithSpan::new(
                SmilesError::UnexpectedLeftParentheses,
                start,
                end,
            ));
        };
        self.push_stack(anchor);
        Ok(())
    }
    /// Evaluates that a right token parentheses has a valid anchor atom.
    ///
    /// # Errors
    /// - Returns [`SmilesError::UnexpectedRightParentheses`] if a valid anchor
    ///   is not found to associate with the right parentheses.
    pub fn validate_branch_close(
        &mut self,
        start: usize,
        end: usize,
    ) -> Result<(), SmilesErrorWithSpan> {
        let Some(anchor) = self.pop_branch_stack() else {
            return Err(SmilesErrorWithSpan::new(
                SmilesError::UnexpectedRightParentheses,
                start,
                end,
            ));
        };
        if self.branch_status != (true, true) {
            return Err(SmilesErrorWithSpan::new(SmilesError::InvalidBranch, start, end));
        }
        self.branch_status = (false, false);
        self.update_last_atom(Some(anchor));
        Ok(())
    }
    /// Checks that there is an existing atom before the current bond, then
    /// updates the pending bond field with the bond.
    ///
    /// # Errors
    /// - Returns [`SmilesError::IncompleteBond`] if a previous atom is not
    ///   found.
    /// - Returns [`SmilesError::InvalidBond`] if bond is not binding two valid
    ///   nodes
    pub fn validate_and_add_bond(
        &mut self,
        start: usize,
        end: usize,
        bond: Bond,
        next_token: Option<&TokenWithSpan>,
    ) -> Result<(), SmilesErrorWithSpan> {
        if self.last_atom().is_none() {
            return Err(SmilesErrorWithSpan::new(SmilesError::IncompleteBond(bond), start, end));
        }
        if let Some(token) = next_token
            && (token.is_bond() || token.token() == Token::LeftParentheses)
        {
            return Err(SmilesErrorWithSpan::new(SmilesError::InvalidBond, start, end));
        }
        self.update_pending_bond(Some(bond));
        Ok(())
    }
    /// Validates that a [`Token::NonBond`] is preceded and proceeded by valid
    /// tokens
    ///
    /// # Errors
    /// - Returns [`SmilesError::InvalidNonBondToken`] if there isn't a valid
    ///   token before or after the non bond
    pub fn validate_non_bond(
        &self,
        last_token: Option<&TokenWithSpan>,
        next_token: Option<&TokenWithSpan>,
        start: usize,
        end: usize,
    ) -> Result<(), SmilesErrorWithSpan> {
        if let Some(last) = last_token {
            match last.token() {
                Token::NonBond
                | Token::BracketedAtom(_)
                | Token::Bond(_)
                | Token::LeftParentheses => {
                    return Err(SmilesErrorWithSpan::new(
                        SmilesError::InvalidNonBondToken,
                        start,
                        end,
                    ));
                }
                _ => {}
            }
        } else {
            return Err(SmilesErrorWithSpan::new(SmilesError::InvalidNonBondToken, start, end));
        }
        if let Some(next) = next_token {
            match next.token() {
                Token::UnbracketedAtom(_) => {}
                _ => {
                    return Err(SmilesErrorWithSpan::new(
                        SmilesError::InvalidNonBondToken,
                        start,
                        end,
                    ));
                }
            }
        } else {
            return Err(SmilesErrorWithSpan::new(SmilesError::InvalidNonBondToken, start, end));
        }
        Ok(())
    }
}

impl Default for ParserState {
    fn default() -> Self {
        Self::new()
    }
}

/// Contains the slice of tokens being iterated on and current position in that
/// slice.
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
        let mut parser_state = ParserState::new();

        while let Some(token_with_span) = self.current() {
            let (start, end) = (token_with_span.start(), token_with_span.end());
            parser_state.update_last_span((start, end));
            match token_with_span.token() {
                Token::UnbracketedAtom(atom) => {
                    parser_state.add_atom(Atom::from(atom), token_with_span.span(), start, end)?;
                }
                Token::BracketedAtom(atom) => {
                    parser_state.add_atom(Atom::from(atom), token_with_span.span(), start, end)?;
                }
                Token::Bond(bond) => {
                    parser_state.validate_and_add_bond(start, end, bond, self.peek_next())?;
                }
                Token::LeftParentheses => {
                    parser_state.validate_branch_open(start, end, self.peek_next())?;
                }
                Token::NonBond => {
                    parser_state.validate_all_closed()?;
                    parser_state.validate_non_bond(
                        self.tokens.get(self.position.saturating_sub(1)),
                        self.tokens.get(self.position.saturating_add(1)),
                        start,
                        end,
                    )?;
                }
                Token::RingClosure(ring_num) => {
                    parser_state.validate_and_add_ring_num(start, end, ring_num)?;
                }
                Token::RightParentheses => {
                    parser_state.validate_branch_close(start, end)?;
                }
            }
            self.advance();
        }
        parser_state.validate_all_closed()?;

        Ok(parser_state.into_smiles())
    }
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
        Smiles, SmilesError,
        atom::{
            Atom, atom_node::AtomNode, atom_symbol::AtomSymbol, bracketed::BracketAtom,
            unbracketed::UnbracketedAtom,
        },
        bond::{Bond, ring_num::RingNum},
        parser::smiles_parser::{ParserState, SmilesParser, default_bond},
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
    }

    #[test]
    fn smiles_parser_advance_moves_position() {
        let tokens = test_tokens();
        let mut p = SmilesParser::new(&tokens);

        p.advance();
        assert_eq!(p.position(), 1);
        assert_eq!(p.current(), Some(&tokens[1]));

        let tokens = vec![TokenWithSpan::new(Token::NonBond, 0, 1)];
        let mut single_token = SmilesParser::new(&tokens);
        single_token.advance();
        single_token.advance();
        assert_eq!(single_token.position(), 1);

        let tokens = vec![TokenWithSpan::new(Token::NonBond, 0, 1)];
        let mut single_token = SmilesParser::new(&tokens);
        single_token.next_token();
        single_token.next_token();
        assert_eq!(single_token.next_token(), None);
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

    #[test]
    fn smiles_parser_errors() {
        use crate::errors::SmilesError;

        fn c(start: usize, end: usize) -> TokenWithSpan {
            TokenWithSpan::new(
                Token::UnbracketedAtom(UnbracketedAtom::new(
                    AtomSymbol::Element(Element::C),
                    false,
                )),
                start,
                end,
            )
        }

        fn ring(start: usize, end: usize, n: u8) -> TokenWithSpan {
            TokenWithSpan::new(Token::RingClosure(RingNum::try_new(n).unwrap()), start, end)
        }

        fn bond(start: usize, end: usize, bond: Bond) -> TokenWithSpan {
            TokenWithSpan::new(Token::Bond(bond), start, end)
        }

        fn lparen(start: usize, end: usize) -> TokenWithSpan {
            TokenWithSpan::new(Token::LeftParentheses, start, end)
        }

        fn rparen(start: usize, end: usize) -> TokenWithSpan {
            TokenWithSpan::new(Token::RightParentheses, start, end)
        }

        fn dot(start: usize, end: usize) -> TokenWithSpan {
            TokenWithSpan::new(Token::NonBond, start, end)
        }

        fn assert_parse_error(
            tokens: &[TokenWithSpan],
            expected_error: SmilesError,
            expected_start: usize,
            expected_end: usize,
        ) {
            let err = SmilesParser::new(tokens).parse().expect_err("expected parser to fail");

            assert_eq!(err.smiles_error(), expected_error);
            assert_eq!(err.start(), expected_start);
            assert_eq!(err.end(), expected_end);
            assert_eq!(err.span(), (expected_start..expected_end));
            assert_eq!(
                err.to_string(),
                format!("{expected_error} at {expected_start}..{expected_end}")
            );
        }

        assert_parse_error(&[lparen(0, 1)], SmilesError::UnexpectedLeftParentheses, 0, 1);

        assert_parse_error(&[rparen(0, 1)], SmilesError::UnexpectedRightParentheses, 0, 1);

        assert_parse_error(&[ring(0, 1, 1)], SmilesError::InvalidRingNumber, 0, 1);

        assert_parse_error(
            &[c(0, 1), bond(1, 2, Bond::Triple), dot(2, 3)],
            SmilesError::IncompleteBond(Bond::Triple),
            2,
            3,
        );

        assert_parse_error(&[c(0, 1), lparen(1, 2), dot(2, 3)], SmilesError::UnclosedBranch, 2, 3);

        assert_parse_error(&[c(0, 1), ring(1, 2, 1), dot(2, 3)], SmilesError::UnclosedRing, 2, 3);

        // parse_end_check branches
        assert_parse_error(
            &[c(0, 1), bond(1, 2, Bond::Double)],
            SmilesError::IncompleteBond(Bond::Double),
            1,
            2,
        );

        assert_parse_error(&[c(0, 1), lparen(1, 2)], SmilesError::UnclosedBranch, 1, 2);

        assert_parse_error(&[c(0, 1), ring(1, 2, 1)], SmilesError::UnclosedRing, 1, 2);
    }

    fn atom_node(
        id: usize,
        element: Element,
        aromatic: bool,
        start: usize,
        end: usize,
    ) -> AtomNode {
        AtomNode::new(
            UnbracketedAtom::new(AtomSymbol::Element(element), aromatic).into(),
            id,
            start..end,
        )
    }

    fn bracket_atom(element: Element) -> BracketAtom {
        BracketAtom::builder().with_symbol(AtomSymbol::Element(element)).build()
    }

    #[test]
    fn parser_state_new_and_default_are_empty() {
        let state = ParserState::new();
        assert_eq!(state.last_span(), (0, 0));
        assert_eq!(state.next_id(), 0);
        assert_eq!(state.last_atom(), None);
        assert_eq!(state.pending_bond(), None);
        assert!(state.branch_stack().is_empty());
        assert!(state.stack_empty());
        assert!(state.ring_open_empty());
        assert!(state.smiles().nodes().is_empty());
        assert!(state.smiles().edges().is_empty());

        let default_state = ParserState::default();
        assert_eq!(default_state.last_span(), (0, 0));
        assert_eq!(default_state.next_id(), 0);
        assert_eq!(default_state.last_atom(), None);
        assert_eq!(default_state.pending_bond(), None);
        assert!(default_state.branch_stack().is_empty());
        assert!(default_state.stack_empty());
        assert!(default_state.ring_open_empty());
        assert!(default_state.smiles().nodes().is_empty());
        assert!(default_state.smiles().edges().is_empty());
    }

    #[test]
    fn parser_state_span_and_id_updates_work() {
        let mut state = ParserState::new();

        state.update_last_span((3, 7));
        assert_eq!(state.last_span(), (3, 7));

        assert_eq!(state.next_id(), 0);
        state.increment_next_id();
        assert_eq!(state.next_id(), 1);
        state.increment_next_id();
        assert_eq!(state.next_id(), 2);
    }

    #[test]
    fn parser_state_last_atom_and_pending_bond_updates_work() {
        let mut state = ParserState::new();

        state.update_last_atom(Some(4));
        assert_eq!(state.last_atom(), Some(4));
        state.update_last_atom(None);
        assert_eq!(state.last_atom(), None);

        state.update_pending_bond(Some(Bond::Triple));
        assert_eq!(state.pending_bond(), Some(Bond::Triple));
        state.update_pending_bond(None);
        assert_eq!(state.pending_bond(), None);
    }

    #[test]
    fn parser_state_branch_stack_methods_work() {
        let mut state = ParserState::new();

        assert!(state.stack_empty());
        assert_eq!(state.branch_stack(), &[]);

        state.push_stack(1);
        state.push_stack(3);

        assert!(!state.stack_empty());
        assert_eq!(state.branch_stack(), &[1, 3]);
        assert_eq!(state.pop_branch_stack(), Some(3));
        assert_eq!(state.pop_branch_stack(), Some(1));
        assert_eq!(state.pop_branch_stack(), None);
        assert!(state.stack_empty());
    }

    #[test]
    fn parser_state_ring_open_methods_work() {
        let mut state = ParserState::new();
        let ring = RingNum::try_new(7).unwrap();

        assert!(state.ring_open_empty());
        assert_eq!(state.remove_ring_open(&ring), None);

        state.insert_ring(ring, (9, Some(Bond::Double)));
        assert!(!state.ring_open_empty());
        assert_eq!(state.remove_ring_open(&ring), Some((9, Some(Bond::Double))));
        assert!(state.ring_open_empty());
    }

    #[test]
    fn parser_state_push_node_and_push_edge_work() {
        let mut state = ParserState::new();

        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();
        state.push_node(atom_node(1, Element::O, false, 1, 2)).unwrap();
        state.push_edge(0, 1, Bond::Double, None).unwrap();

        assert_eq!(state.smiles().nodes().len(), 2);
        assert_eq!(state.smiles().edges().len(), 1);
        assert!(state.edge_for_node_pair_exists((0, 1)));
    }

    #[test]
    fn parser_state_push_node_duplicate_errors() {
        let mut state = ParserState::new();

        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();
        let err = state
            .push_node(atom_node(0, Element::O, false, 1, 2))
            .expect_err("expected duplicate node id");

        assert_eq!(err, SmilesError::DuplicateNodeId(0));
    }

    #[test]
    fn parser_state_push_edge_invalid_node_errors() {
        let mut state = ParserState::new();
        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();

        let err_a =
            state.push_edge(9, 0, Bond::Single, None).expect_err("expected invalid first node");
        assert_eq!(err_a, SmilesError::NodeIdInvalid(9));

        let err_b =
            state.push_edge(0, 8, Bond::Single, None).expect_err("expected invalid second node");
        assert_eq!(err_b, SmilesError::NodeIdInvalid(8));
    }

    #[test]
    fn parser_state_into_smiles_returns_owned_graph() {
        let mut state = ParserState::new();
        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();
        state.push_node(atom_node(1, Element::O, false, 1, 2)).unwrap();
        state.push_edge(0, 1, Bond::Single, None).unwrap();

        let smiles = state.into_smiles();
        assert_eq!(smiles.nodes().len(), 2);
        assert_eq!(smiles.edges().len(), 1);
    }

    #[test]
    fn parser_state_add_atom_adds_first_atom() {
        let mut state = ParserState::new();

        state
            .add_atom(
                Atom::from(UnbracketedAtom::new(AtomSymbol::Element(Element::C), false)),
                0..1,
                0,
                1,
            )
            .unwrap();

        assert_eq!(state.next_id(), 1);
        assert_eq!(state.last_atom(), Some(0));
        assert_eq!(state.pending_bond(), None);
        assert_eq!(state.smiles().nodes().len(), 1);
        assert_eq!(state.smiles().edges().len(), 0);
    }

    #[test]
    fn parser_state_add_atom_connects_with_pending_bond() {
        let mut state = ParserState::new();

        state
            .add_atom(
                Atom::from(UnbracketedAtom::new(AtomSymbol::Element(Element::C), false)),
                0..1,
                0,
                1,
            )
            .unwrap();
        state.update_pending_bond(Some(Bond::Triple));
        state
            .add_atom(
                Atom::from(UnbracketedAtom::new(AtomSymbol::Element(Element::N), false)),
                1..2,
                1,
                2,
            )
            .unwrap();

        assert_eq!(state.next_id(), 2);
        assert_eq!(state.last_atom(), Some(1));
        assert_eq!(state.pending_bond(), None);
        assert_eq!(state.smiles().nodes().len(), 2);
        assert_eq!(state.smiles().edges().len(), 1);
        assert_eq!(state.smiles().edges()[0].bond(), Bond::Triple);
    }

    #[test]
    fn parser_state_validate_all_closed_ok_resets_last_atom_and_pending_bond() {
        let mut state = ParserState::new();
        state.update_last_span((4, 5));
        state.update_last_atom(Some(9));
        state.update_pending_bond(None);

        state.validate_all_closed().unwrap();

        assert_eq!(state.last_atom(), None);
        assert_eq!(state.pending_bond(), None);
    }

    #[test]
    fn parser_state_validate_all_closed_errors_for_incomplete_bond() {
        let mut state = ParserState::new();
        state.update_last_span((2, 3));
        state.update_pending_bond(Some(Bond::Double));

        let err = state.validate_all_closed().expect_err("expected incomplete bond");

        assert_eq!(err.smiles_error(), SmilesError::IncompleteBond(Bond::Double));
        assert_eq!(err.start(), 2);
        assert_eq!(err.end(), 3);
    }

    #[test]
    fn parser_state_validate_all_closed_errors_for_unclosed_branch() {
        let mut state = ParserState::new();
        state.update_last_span((2, 3));
        state.push_stack(0);

        let err = state.validate_all_closed().expect_err("expected unclosed branch");

        assert_eq!(err.smiles_error(), SmilesError::UnclosedBranch);
        assert_eq!(err.start(), 2);
        assert_eq!(err.end(), 3);
    }

    #[test]
    fn parser_state_validate_all_closed_errors_for_unclosed_ring() {
        let mut state = ParserState::new();
        state.update_last_span((2, 3));
        state.insert_ring(RingNum::try_new(1).unwrap(), (0, None));

        let err = state.validate_all_closed().expect_err("expected unclosed ring");

        assert_eq!(err.smiles_error(), SmilesError::UnclosedRing);
        assert_eq!(err.start(), 2);
        assert_eq!(err.end(), 3);
    }

    #[test]
    fn parser_state_validate_branch_open_errors_without_anchor() {
        let mut state = ParserState::new();

        let err = state.validate_branch_open(3, 4, None).expect_err("expected missing anchor");

        assert_eq!(err.smiles_error(), SmilesError::UnexpectedLeftParentheses);
        assert_eq!(err.start(), 3);
        assert_eq!(err.end(), 4);
    }

    #[test]
    fn parser_state_validate_branch_close_errors_without_open_branch() {
        let mut state = ParserState::new();

        let err = state.validate_branch_close(4, 5).expect_err("expected missing branch anchor");

        assert_eq!(err.smiles_error(), SmilesError::UnexpectedRightParentheses);
        assert_eq!(err.start(), 4);
        assert_eq!(err.end(), 5);
    }

    #[test]
    fn parser_state_validate_and_add_bond_sets_pending_bond() {
        let mut state = ParserState::new();
        state.update_last_atom(Some(0));

        state.validate_and_add_bond(1, 2, Bond::Aromatic, None).unwrap();

        assert_eq!(state.pending_bond(), Some(Bond::Aromatic));
    }

    #[test]
    fn parser_state_validate_and_add_bond_errors_without_left_atom() {
        let mut state = ParserState::new();

        let err = state
            .validate_and_add_bond(1, 2, Bond::Single, None)
            .expect_err("expected incomplete bond");

        assert_eq!(err.smiles_error(), SmilesError::IncompleteBond(Bond::Single));
        assert_eq!(err.start(), 1);
        assert_eq!(err.end(), 2);
    }

    #[test]
    fn parser_state_validate_and_add_ring_num_opens_ring_and_clears_pending_bond() {
        let mut state = ParserState::new();
        let ring = RingNum::try_new(6).unwrap();

        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();
        state.update_last_atom(Some(0));
        state.update_pending_bond(Some(Bond::Double));

        state.validate_and_add_ring_num(1, 2, ring).unwrap();

        assert_eq!(state.pending_bond(), None);
        assert_eq!(state.remove_ring_open(&ring), Some((0, Some(Bond::Double))));
    }

    #[test]
    fn parser_state_validate_and_add_ring_num_closes_ring_and_adds_edge() {
        let mut state = ParserState::new();
        let ring = RingNum::try_new(4).unwrap();

        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();
        state.push_node(atom_node(1, Element::N, false, 1, 2)).unwrap();
        state.insert_ring(ring, (0, Some(Bond::Triple)));
        state.update_last_atom(Some(1));

        state.validate_and_add_ring_num(2, 3, ring).unwrap();

        assert_eq!(state.pending_bond(), None);
        assert!(state.edge_for_node_pair_exists((0, 1)));
        assert_eq!(state.smiles().edges().len(), 1);
        assert_eq!(state.smiles().edges()[0].bond(), Bond::Triple);
        assert_eq!(state.smiles().edges()[0].ring_num_val(), Some(4));
    }

    #[test]
    fn parser_state_validate_and_add_ring_num_prefers_current_pending_bond() {
        let mut state = ParserState::new();
        let ring = RingNum::try_new(5).unwrap();

        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();
        state.push_node(atom_node(1, Element::O, false, 1, 2)).unwrap();
        state.insert_ring(ring, (0, Some(Bond::Double)));
        state.update_last_atom(Some(1));
        state.update_pending_bond(Some(Bond::Quadruple));

        state.validate_and_add_ring_num(2, 3, ring).unwrap();

        assert_eq!(state.smiles().edges()[0].bond(), Bond::Quadruple);
    }

    #[test]
    fn parser_state_validate_and_add_ring_num_errors_without_current_atom() {
        let mut state = ParserState::new();
        let ring = RingNum::try_new(1).unwrap();

        let err =
            state.validate_and_add_ring_num(7, 8, ring).expect_err("expected invalid ring number");

        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
        assert_eq!(err.start(), 7);
        assert_eq!(err.end(), 8);
    }

    #[test]
    fn parser_state_validate_and_add_ring_num_errors_for_self_loop() {
        let mut state = ParserState::new();
        let ring = RingNum::try_new(2).unwrap();

        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();
        state.insert_ring(ring, (0, None));
        state.update_last_atom(Some(0));

        let err = state
            .validate_and_add_ring_num(3, 4, ring)
            .expect_err("expected invalid self-loop ring");

        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
        assert_eq!(err.start(), 3);
        assert_eq!(err.end(), 4);
    }

    #[test]
    fn parser_state_validate_and_add_ring_num_errors_for_duplicate_edge() {
        let mut state = ParserState::new();
        let ring = RingNum::try_new(3).unwrap();

        state.push_node(atom_node(0, Element::C, false, 0, 1)).unwrap();
        state.push_node(atom_node(1, Element::O, false, 1, 2)).unwrap();
        state.push_edge(0, 1, Bond::Single, None).unwrap();
        state.insert_ring(ring, (0, None));
        state.update_last_atom(Some(1));

        let err = state
            .validate_and_add_ring_num(4, 5, ring)
            .expect_err("expected duplicate edge ring failure");

        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
        assert_eq!(err.start(), 4);
        assert_eq!(err.end(), 5);
    }

    #[test]
    fn default_bond_is_aromatic_for_two_aromatic_atoms() {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C, true, 0, 1)).unwrap();
        smiles.push_node(atom_node(1, Element::N, true, 1, 2)).unwrap();

        assert_eq!(default_bond(&smiles, 0, 1), Bond::Aromatic);
    }

    #[test]
    fn default_bond_is_single_when_atoms_are_not_both_aromatic() {
        let mut smiles = Smiles::new();
        smiles.push_node(atom_node(0, Element::C, true, 0, 1)).unwrap();
        smiles.push_node(atom_node(1, Element::O, false, 1, 2)).unwrap();

        assert_eq!(default_bond(&smiles, 0, 1), Bond::Single);
    }

    #[test]
    fn smiles_parser_parse_builds_graph_for_bracketed_atom() {
        let tokens = vec![TokenWithSpan::new(Token::BracketedAtom(bracket_atom(Element::N)), 0, 3)];

        let smiles = SmilesParser::new(&tokens).parse().expect("parse should succeed");

        assert_eq!(smiles.nodes().len(), 1);
        assert_eq!(smiles.edges().len(), 0);
    }

    #[test]
    fn smiles_parser_rejects_self_loop_ring_closure() {
        let tokens = vec![
            TokenWithSpan::new(
                Token::UnbracketedAtom(UnbracketedAtom::new(
                    AtomSymbol::Element(Element::C),
                    false,
                )),
                0,
                1,
            ),
            TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1).unwrap()), 1, 2),
            TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1).unwrap()), 2, 3),
        ];

        let err = SmilesParser::new(&tokens).parse().expect_err("expected self-loop ring failure");

        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
        assert_eq!(err.start(), 2);
        assert_eq!(err.end(), 3);
    }

    #[test]
    fn smiles_parser_rejects_duplicate_ring_edge_between_same_atoms() {
        let tokens = vec![
            TokenWithSpan::new(
                Token::UnbracketedAtom(UnbracketedAtom::new(
                    AtomSymbol::Element(Element::C),
                    false,
                )),
                0,
                1,
            ),
            TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1).unwrap()), 1, 2),
            TokenWithSpan::new(
                Token::UnbracketedAtom(UnbracketedAtom::new(
                    AtomSymbol::Element(Element::O),
                    false,
                )),
                2,
                3,
            ),
            TokenWithSpan::new(Token::RingClosure(RingNum::try_new(1).unwrap()), 3, 4),
        ];

        let err =
            SmilesParser::new(&tokens).parse().expect_err("expected duplicate edge ring failure");

        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
        assert_eq!(err.start(), 3);
        assert_eq!(err.end(), 4);
    }
}
