//! Second pass that parses the [`TokenWithSpan`]

use std::{collections::HashMap, ops::Range};

use crate::{
    atom::{Atom, atom_node::AtomNode},
    bond::{Bond, ring_num::RingNum},
    errors::{SmilesError, SmilesErrorWithSpan},
    smiles::Smiles,
    token::{Token, TokenWithSpan},
};

/// Structure for containing the Parser State
pub struct ParserState {
    smiles: Smiles,
    next_id: usize,
    last_atom: Option<usize>,
    pending_bond: Option<Bond>,
    branch_stack: Vec<usize>,
    ring_open: HashMap<RingNum, (usize, Option<Bond>)>,
    last_span: (usize, usize),
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
    /// Returns the next id field.
    #[must_use]
    pub fn next_id(&self) -> usize {
        self.next_id
    }
    /// Updates the next available id by incrementing.
    pub fn increment_next_id(&mut self) {
        self.next_id += 1;
    }
    /// Returns the previous atom id if there is one.
    #[must_use]
    pub fn last_atom(&self) -> Option<usize> {
        self.last_atom
    }
    /// Updates the last item id.
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
    /// Iserts the give ring ito the ring open field
    pub fn insert_ring(&mut self, ring_num: RingNum, pending: (usize, Option<Bond>)) {
        self.ring_open.insert(ring_num, pending);
    }
    /// Returns a borrowed reference of the smiles field.
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
    /// - Returns [`SmilesError::DuplicateNodeId`] if node id already exists
    pub fn push_node(&mut self, node: AtomNode) -> Result<(), SmilesError> {
        self.smiles.push_node(node)
    }
    /// Pushes an [`BondEdge`] into the parsed [`Smiles`] graph.
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
    /// Adds an atom to the smiles graph, either bracketed or unbracketed. 
    pub fn add_atom(&mut self, atom: Atom, token_span: Range<usize>, start: usize, end: usize) -> Result<(), SmilesErrorWithSpan>{
        let id = self.next_id();
        self.increment_next_id();
        let node = AtomNode::new(atom, id, token_span);
        self.push_node(node).map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;
        if let Some(prev) = self.last_atom() {
            let bond = self.pending_bond().unwrap_or_else(|| default_bond(self.smiles(), prev, id));
            self.push_edge(prev, id, bond, None).map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;
        }
        self.update_last_atom(Some(id));
        self.update_pending_bond(None);
        Ok(())
    }
}

impl Default for ParserState {
    fn default() -> Self {
        Self::new()
    }
}

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
                    parser_state.add_atom(Atom::Unbracketed(atom), token_with_span.span(), start, end)?;
                }
                Token::BracketedAtom(atom) => {
                    parser_state.add_atom(Atom::Bracketed(atom), token_with_span.span(), start, end)?;
                }
                Token::Bond(bond) => {
                    if parser_state.last_atom().is_none() {
                        return Err(SmilesErrorWithSpan::new(
                            SmilesError::IncompleteBond(bond),
                            start,
                            end,
                        ));
                    }
                    parser_state.update_pending_bond(Some(bond));
                }
                Token::LeftParentheses => {
                    let Some(anchor) = parser_state.last_atom() else {
                        return Err(SmilesErrorWithSpan::new(
                            SmilesError::UnexpectedLeftParentheses,
                            start,
                            end,
                        ));
                    };
                    parser_state.push_stack(anchor);
                }
                Token::NonBond => {
                    validate_non_bond(
                        parser_state.pending_bond(),
                        parser_state.branch_stack().is_empty(),
                        parser_state.ring_open_empty(),
                        parser_state.last_span(),
                    )?;
                    parser_state.update_last_atom(None);
                    parser_state.update_pending_bond(None);
                }
                Token::RingClosure(ring_num) => {
                    let Some(current) = parser_state.last_atom() else {
                        return Err(SmilesErrorWithSpan::new(
                            SmilesError::InvalidRingNumber,
                            start,
                            end,
                        ));
                    };
                    if let Some((other, stored_bond)) = parser_state.remove_ring_open(&ring_num) {
                        if current == other {
                            return Err(SmilesErrorWithSpan::new(
                                SmilesError::InvalidRingNumber,
                                start,
                                end,
                            ));
                        }
                        if parser_state.edge_for_node_pair_exists((current, other)) {
                            return Err(SmilesErrorWithSpan::new(
                                SmilesError::InvalidRingNumber,
                                start,
                                end,
                            ));
                        }
                        let bond = parser_state
                            .pending_bond()
                            .or(stored_bond)
                            .unwrap_or_else(|| default_bond(parser_state.smiles(), current, other));

                        parser_state
                            .push_edge(current, other, bond, Some(ring_num))
                            .map_err(|e| SmilesErrorWithSpan::new(e, start, end))?;

                        parser_state.update_pending_bond(None);
                    } else {
                        parser_state.insert_ring(ring_num, (current, parser_state.pending_bond));
                        parser_state.update_pending_bond(None);
                    }
                }
                Token::RightParentheses => {
                    let Some(anchor) = parser_state.pop_branch_stack() else {
                        return Err(SmilesErrorWithSpan::new(
                            SmilesError::UnexpectedRightParentheses,
                            start,
                            end,
                        ));
                    };
                    parser_state.update_last_atom(Some(anchor));
                }
            }
            self.advance();
        }
        parse_end_check(
            parser_state.pending_bond(),
            parser_state.stack_empty(),
            parser_state.ring_open_empty(),
            parser_state.last_span(),
        )?;

        Ok(parser_state.into_smiles())
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
}
