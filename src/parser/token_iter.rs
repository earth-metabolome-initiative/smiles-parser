//! Submodule creating the `TokenIter` struct, which is an iterator over
//! the `Token`s found in a provided string.

use std::str::FromStr;

use elements_rs::Element;

use crate::{
    atom::{
        atom_symbol::AtomSymbol,
        bracketed::{
            BracketAtom, charge::Charge, chirality::Chirality, hydrogen_count::HydrogenCount,
        },
        unbracketed::UnbracketedAtom,
    },
    bond::{Bond, ring_num::RingNum},
    errors::{SmilesError, SmilesErrorWithSpan},
    token::{Token, TokenWithSpan},
};

/// An iterator over the tokens found in a SMILES string.
pub struct TokenIter<'a> {
    /// The peekable `Chars` with `Indices` iterator
    chars: std::iter::Peekable<std::str::CharIndices<'a>>,
    /// Denotes whether currently inside brackets
    in_bracket: bool,
    /// The length of the input
    len: usize,
}

impl<'a> From<&'a str> for TokenIter<'a> {
    fn from(s: &'a str) -> Self {
        TokenIter { chars: s.char_indices().peekable(), in_bracket: false, len: s.len() }
    }
}

impl TokenIter<'_> {
    fn parse_token(&mut self, current_char: char) -> Result<Token, SmilesError> {
        let token = match current_char {
            '.' => {
                if self.in_bracket {
                    return Err(SmilesError::NonBondInBracket);
                }
                Token::NonBond
            }
            '[' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedLeftBracket);
                }
                self.in_bracket = true;
                let mut possible_bracket_atom = BracketAtom::builder();
                if let Some(isotope) = try_fold_number(self) {
                    possible_bracket_atom = possible_bracket_atom.with_isotope(isotope?);
                }
                let (atom, aromatic) = try_element(self)?;
                possible_bracket_atom =
                    possible_bracket_atom.with_symbol(atom).with_aromatic(aromatic);
                if let Some(chiral) = try_chirality(self)? {
                    possible_bracket_atom = possible_bracket_atom.with_chiral(chiral);
                }

                if possible_bracket_atom.symbol() == AtomSymbol::Unspecified {
                    return Err(SmilesError::MissingBracketElement);
                }
                possible_bracket_atom = possible_bracket_atom.with_hydrogens(hydrogen_count(self)?);
                possible_bracket_atom = possible_bracket_atom.with_charge(try_charge(self)?);
                possible_bracket_atom = possible_bracket_atom.with_class(try_class(self)?);
                let bracket_atom = possible_bracket_atom.build();
                if matches!(self.peek_char(), Some(']')) {
                    self.in_bracket = false;
                    self.chars.next();
                    Token::BracketedAtom(bracket_atom)
                } else {
                    return Err(SmilesError::UnclosedBracket);
                }
            }
            c if c.is_ascii_alphabetic() || c == '*' => {
                let (symbol, aromatic) = try_element_from_first(self, c)?;
                if !valid_unbracketed(symbol) {
                    return Err(SmilesError::InvalidUnbracketedAtom(symbol));
                }
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedBracketedState);
                }
                Token::UnbracketedAtom(UnbracketedAtom::new(symbol, aromatic))
            }

            n if n.is_ascii_digit() || n == '%' => {
                if n == '%' {
                    if self.in_bracket {
                        return Err(SmilesError::UnexpectedPercent);
                    }

                    if let Some(num) = try_fold_number::<u8>(self) {
                        let ring_num = RingNum::try_new(num?)?;
                        if ring_num.get() < 10 {
                            return Err(SmilesError::InvalidRingNumber);
                        }
                        Token::RingClosure(ring_num)
                    } else {
                        return Err(SmilesError::InvalidRingNumber);
                    }
                } else {
                    let Some(first) = n.to_digit(10) else {
                        return Err(SmilesError::InvalidClass);
                    };

                    Token::RingClosure(RingNum::try_new(u8::try_from(first)?)?)
                }
            }
            '-' | '=' | '#' | '$' | ':' | '/' | '\\' => try_bond(current_char, self.in_bracket)?,
            '(' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedBracketedState);
                }
                Token::LeftParentheses
            }
            ')' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedBracketedState);
                }
                Token::RightParentheses
            }
            _ => return Err(SmilesError::UnexpectedCharacter(current_char)),
        };
        Ok(token)
    }

    fn current_end(&mut self) -> usize {
        if let Some(&(next_id, _)) = self.chars.peek() { next_id } else { self.len }
    }
    fn peek_char(&mut self) -> Option<char> {
        self.chars.peek().map(|(_, c)| *c)
    }
    fn next_char(&mut self) -> Option<char> {
        self.chars.next().map(|(_, c)| c)
    }
}

impl Iterator for TokenIter<'_> {
    type Item = Result<TokenWithSpan, SmilesErrorWithSpan>;

    fn next(&mut self) -> Option<Self::Item> {
        let (start, current_char) = self.chars.next()?;
        match self.parse_token(current_char) {
            Ok(token) => {
                let end = self.current_end();
                Some(Ok(TokenWithSpan::new(token, start, end)))
            }
            Err(e) => {
                let mut end = self.current_end();
                if end <= start {
                    end = (start + current_char.len_utf8()).min(self.len);
                }
                Some(Err(SmilesErrorWithSpan::new(e, start, end)))
            }
        }
    }
}

/// determines whether an aromatic is valid for a given bracketed or unbracketed
/// atom
///
/// # Parameters
/// - `bool` for the status of `in_bracket`
/// - the [`Element`] being passed
fn aromatic_from_element(in_bracket: bool, element: Element) -> Result<bool, SmilesError> {
    let allowed = if in_bracket {
        matches!(
            element,
            Element::B
                | Element::C
                | Element::N
                | Element::O
                | Element::P
                | Element::S
                | Element::Se
                | Element::As
        )
    } else {
        matches!(
            element,
            Element::B | Element::C | Element::N | Element::O | Element::S | Element::P
        )
    };
    if allowed { Ok(true) } else { Err(SmilesError::InvalidAromaticElement(element)) }
}

fn try_element(stream: &mut TokenIter<'_>) -> Result<(AtomSymbol, bool), SmilesError> {
    let first = stream.next_char().ok_or(SmilesError::MissingElement)?;
    try_element_from_first(stream, first)
}

fn try_element_from_first(
    stream: &mut TokenIter<'_>,
    char_1: char,
) -> Result<(AtomSymbol, bool), SmilesError> {
    if char_1 == '*' {
        return Ok((AtomSymbol::WildCard, false));
    }
    if !char_1.is_ascii_alphabetic() {
        return Err(SmilesError::MissingElement);
    }

    let is_aromatic_candidate = char_1.is_ascii_lowercase();
    let try_candidate = |val: &str| -> Option<Element> { Element::from_str(val).ok() };

    if let Some(char_2) = stream.peek_char()
        && char_2.is_ascii_alphabetic()
    {
        if is_aromatic_candidate && char_2.is_ascii_lowercase() {
            let candidate = format!("{}{}", char_1.to_ascii_uppercase(), char_2);
            if let Some(element) = try_candidate(&candidate) {
                stream.chars.next();
                let aromatic = aromatic_from_element(stream.in_bracket, element)?;
                return Ok((AtomSymbol::Element(element), aromatic));
            }
        }
        if !is_aromatic_candidate && char_2.is_ascii_lowercase() {
            let candidate = format!("{char_1}{char_2}");
            if let Some(element) = try_candidate(&candidate) {
                stream.chars.next();
                return Ok((AtomSymbol::Element(element), false));
            }
        }
    }

    let one = if is_aromatic_candidate {
        char_1.to_ascii_uppercase().to_string()
    } else {
        char_1.to_string()
    };

    if let Some(element) = try_candidate(&one) {
        let aromatic = if is_aromatic_candidate {
            aromatic_from_element(stream.in_bracket, element)?
        } else {
            false
        };
        return Ok((AtomSymbol::Element(element), aromatic));
    }

    Err(SmilesError::InvalidElementName(char_1))
}

// B, C, N, O, P, S, F, Cl, Br, I,
fn valid_unbracketed(symbol: AtomSymbol) -> bool {
    match symbol {
        AtomSymbol::Element(element) => {
            matches!(
                element,
                Element::B
                    | Element::C
                    | Element::N
                    | Element::O
                    | Element::P
                    | Element::S
                    | Element::F
                    | Element::Cl
                    | Element::Br
                    | Element::I
            )
        }
        AtomSymbol::WildCard => true,
        AtomSymbol::Unspecified => false,
    }
}

fn try_chirality(stream: &mut TokenIter<'_>) -> Result<Option<Chirality>, SmilesError> {
    if stream.peek_char() != Some('@') {
        return Ok(None);
    }
    stream.chars.next();
    let char_2 = stream.peek_char().ok_or(SmilesError::UnexpectedEndOfString)?;
    let chirality = match char_2 {
        '@' => {
            stream.chars.next();
            Chirality::AtAt
        }
        'T' => {
            match stream.peek_char().ok_or(SmilesError::UnexpectedEndOfString)? {
                'H' => {
                    stream.chars.next();
                    let num =
                        try_fold_number::<u8>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_th(num)?
                }
                'B' => {
                    stream.chars.next();
                    let num =
                        try_fold_number::<u8>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_tb(num)?
                }
                _ => return Err(SmilesError::InvalidChirality),
            }
        }
        'A' | 'S' => {
            stream.chars.next();
            match stream.peek_char().ok_or(SmilesError::UnexpectedEndOfString)? {
                'P' => {
                    stream.chars.next();
                    let num =
                        try_fold_number::<u8>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_sp(num)?
                }
                _ => return Err(SmilesError::InvalidChirality),
            }
        }
        'O' => {
            stream.chars.next();
            match stream.peek_char().ok_or(SmilesError::UnexpectedEndOfString)? {
                'H' => {
                    stream.chars.next();
                    let num =
                        try_fold_number::<u8>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_oh(num)?
                }
                _ => return Err(SmilesError::InvalidChirality),
            }
        }
        'H' | '-' | '+' | ':' | ']' => Chirality::At,
        _ => return Err(SmilesError::InvalidChirality),
    };
    Ok(Some(chirality))
}

fn try_fold_number<B>(stream: &mut TokenIter<'_>) -> Option<Result<B, SmilesError>>
where
    B: TryFrom<u32>,
{
    let mut seen_any = false;
    let mut amount: u32 = 0;

    while let Some(char) = stream.peek_char() {
        let Some(digit) = char.to_digit(10) else {
            break;
        };
        stream.chars.next();
        seen_any = true;
        match amount.checked_mul(10).and_then(|x| x.checked_add(digit)) {
            Some(val) => amount = val,
            None => return Some(Err(SmilesError::IntegerOverflow)),
        }
    }

    if !seen_any {
        return None;
    }

    Some(B::try_from(amount).map_err(|_| SmilesError::IntegerOverflow))
}

fn hydrogen_count(stream: &mut TokenIter<'_>) -> Result<HydrogenCount, SmilesError> {
    let possible_hydrogen = stream.peek_char();
    if matches!(possible_hydrogen, Some('H')) {
        stream.chars.next();
        match try_fold_number::<u8>(stream) {
            Some(h) => Ok(HydrogenCount::new(Some(h?))),
            None => Ok(HydrogenCount::new(Some(1))),
        }
    } else {
        Ok(HydrogenCount::Unspecified)
    }
}

fn try_charge(stream: &mut TokenIter<'_>) -> Result<Charge, SmilesError> {
    match stream.peek_char() {
        Some('-') => {
            stream.chars.next();
            match stream.peek_char() {
                Some('-') => {
                    stream.chars.next();
                    Charge::try_new(-2)
                }
                _ => {
                    if let Some(possible_num) = try_fold_number::<i8>(stream) {
                        Charge::try_new(-possible_num?)
                    } else {
                        Charge::try_new(-1)
                    }
                }
            }
        }
        Some('+') => {
            stream.chars.next();
            match stream.peek_char() {
                Some('+') => {
                    stream.chars.next();
                    Charge::try_new(2)
                }
                _ => {
                    if let Some(possible_num) = try_fold_number::<i8>(stream) {
                        Charge::try_new(possible_num?)
                    } else {
                        Charge::try_new(1)
                    }
                }
            }
        }
        _ => Ok(Charge::default()),
    }
}

fn try_class(stream: &mut TokenIter<'_>) -> Result<u16, SmilesError> {
    match stream.peek_char() {
        Some(':') => {
            stream.chars.next();
            if let Some(possible_num) = try_fold_number(stream) {
                possible_num
            } else {
                Err(SmilesError::InvalidClass)
            }
        }
        _ => Ok(0),
    }
}

fn try_bond(char: char, bracket: bool) -> Result<Token, SmilesError> {
    let bond = match char {
        '-' => {
            if bracket {
                return Err(SmilesError::UnexpectedDash);
            }
            Token::Bond(Bond::Single)
        }
        '=' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Double));
            }
            Token::Bond(Bond::Double)
        }
        '#' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Triple));
            }
            Token::Bond(Bond::Triple)
        }
        '$' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Quadruple));
            }
            Token::Bond(Bond::Quadruple)
        }
        ':' => {
            if bracket {
                return Err(SmilesError::UnexpectedColon);
            }
            Token::Bond(Bond::Aromatic)
        }
        '/' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Up));
            }
            Token::Bond(Bond::Up)
        }
        '\\' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Down));
            }
            Token::Bond(Bond::Down)
        }
        _ => return Err(SmilesError::UnexpectedCharacter(char)),
    };
    Ok(bond)
}
