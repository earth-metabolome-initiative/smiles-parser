//! Submodule creating the `TokenIter` struct, which is an iterator over
//! the `Token`s found in a provided string.

use core::str::{FromStr, from_utf8};

use elements_rs::Element;

use crate::{
    atom::{
        Atom,
        atom_symbol::AtomSymbol,
        bracketed::{charge::Charge, chirality::Chirality},
    },
    bond::{Bond, ring_num::RingNum},
    errors::{SmilesError, SmilesErrorWithSpan},
    token::{Token, TokenWithSpan},
};

/// An iterator over the tokens found in a SMILES string.
pub(crate) struct TokenIter<'a> {
    /// Raw input bytes for the ASCII-heavy parsing fast path.
    bytes: &'a [u8],
    /// Current byte offset in the input.
    position: usize,
    /// Denotes whether currently inside brackets
    in_bracket: bool,
    /// The length of the input
    len: usize,
}

impl<'a> From<&'a str> for TokenIter<'a> {
    #[inline]
    fn from(s: &'a str) -> Self {
        TokenIter { bytes: s.as_bytes(), position: 0, in_bracket: false, len: s.len() }
    }
}

impl TokenIter<'_> {
    #[inline]
    fn parse_token(&mut self, current_byte: u8) -> Result<Token, SmilesError> {
        let token = match current_byte {
            b'.' => {
                if self.in_bracket {
                    return Err(SmilesError::NonBondInBracket);
                }
                Token::NonBond
            }
            b'[' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedLeftBracket);
                }
                self.in_bracket = true;
                let isotope_mass_number = if let Some(isotope) = try_fold_number::<u16, 3>(self) {
                    Some(isotope?)
                } else {
                    None
                };
                let (symbol, aromatic) = try_element(self)?;
                let chirality = try_chirality(self)?;
                let hydrogens = hydrogen_count(self)?;
                let charge = try_charge(self)?;
                let class = try_class(self)?;
                let atom = Atom::new_bracket(
                    symbol,
                    isotope_mass_number,
                    aromatic,
                    hydrogens,
                    charge,
                    class,
                    chirality,
                );
                if self.peek_byte() == Some(b']') {
                    self.in_bracket = false;
                    let _ = self.next_byte();
                    Token::Atom(atom)
                } else {
                    return Err(SmilesError::UnclosedBracket);
                }
            }
            c if c.is_ascii_alphabetic() || c == b'*' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedBracketedState);
                }
                let (symbol, aromatic) = if let Some(atom) = try_organic_subset_from_first(self, c)
                {
                    atom?
                } else {
                    let (symbol, aromatic) = try_element_from_first(self, c)?;
                    if !valid_unbracketed(symbol) {
                        return Err(SmilesError::InvalidUnbracketedAtom(symbol));
                    }
                    (symbol, aromatic)
                };
                Token::Atom(Atom::new_organic_subset(symbol, aromatic))
            }

            n if n.is_ascii_digit() || n == b'%' => {
                if n == b'%' {
                    if self.in_bracket {
                        return Err(SmilesError::UnexpectedPercent);
                    }

                    if let Some(num) = try_fold_number::<u8, 2>(self) {
                        let ring_num = RingNum::try_new(num?)?;
                        if ring_num.get() < 10 {
                            return Err(SmilesError::InvalidRingNumber);
                        }
                        Token::RingClosure(ring_num)
                    } else {
                        return Err(SmilesError::InvalidRingNumber);
                    }
                } else {
                    Token::RingClosure(RingNum::try_new(n - b'0')?)
                }
            }
            b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\' => {
                try_bond(current_byte, self.in_bracket)?
            }
            b'(' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedBracketedState);
                }
                Token::LeftParentheses
            }
            b')' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedBracketedState);
                }
                Token::RightParentheses
            }
            _ => return Err(SmilesError::UnexpectedCharacter(char::from(current_byte))),
        };
        Ok(token)
    }

    #[inline]
    fn current_end(&self) -> usize {
        self.position
    }

    #[inline]
    fn peek_byte(&self) -> Option<u8> {
        self.bytes.get(self.position).copied()
    }

    #[inline]
    fn next_byte(&mut self) -> Option<u8> {
        let byte = self.bytes.get(self.position).copied()?;
        self.position += 1;
        Some(byte)
    }
}

impl Iterator for TokenIter<'_> {
    type Item = Result<TokenWithSpan, SmilesErrorWithSpan>;

    fn next(&mut self) -> Option<Self::Item> {
        let start = self.position;
        let current_byte = self.next_byte()?;
        if !current_byte.is_ascii() {
            self.position = (start + utf8_char_width(current_byte)).min(self.len);
            return Some(Err(SmilesErrorWithSpan::new(
                SmilesError::UnexpectedUnicodeCharacter,
                start,
                self.position,
            )));
        }
        match self.parse_token(current_byte) {
            Ok(token) => {
                let end = self.current_end();
                Some(Ok(TokenWithSpan::new(token, start, end)))
            }
            Err(e) => Some(Err(SmilesErrorWithSpan::new(e, start, self.current_end()))),
        }
    }
}

#[inline]
const fn utf8_char_width(first_byte: u8) -> usize {
    match first_byte {
        0xc0..=0xdf => 2,
        0xe0..=0xef => 3,
        0xf0..=0xf7 => 4,
        _ => 1,
    }
}

/// determines whether an aromatic is valid for a given bracketed or unbracketed
/// atom
///
/// # Parameters
/// - `bool` for the status of `in_bracket`
/// - the [`Element`] being passed
#[inline]
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

#[inline]
fn try_element(stream: &mut TokenIter<'_>) -> Result<(AtomSymbol, bool), SmilesError> {
    let first = stream.next_byte().ok_or(SmilesError::MissingElement)?;
    try_element_from_first(stream, first)
}

#[inline]
fn try_organic_subset_from_first(
    stream: &mut TokenIter<'_>,
    byte_1: u8,
) -> Option<Result<(AtomSymbol, bool), SmilesError>> {
    let element = match byte_1 {
        b'*' => return Some(Ok((AtomSymbol::WildCard, false))),
        b'B' => {
            if stream.peek_byte() == Some(b'r') {
                let _ = stream.next_byte();
                Element::Br
            } else {
                Element::B
            }
        }
        b'C' => {
            if stream.peek_byte() == Some(b'l') {
                let _ = stream.next_byte();
                Element::Cl
            } else {
                Element::C
            }
        }
        b'N' => Element::N,
        b'O' => Element::O,
        b'P' => Element::P,
        b'S' => Element::S,
        b'F' => Element::F,
        b'I' => Element::I,
        b'b' => return Some(Ok((AtomSymbol::Element(Element::B), true))),
        b'c' => return Some(Ok((AtomSymbol::Element(Element::C), true))),
        b'n' => return Some(Ok((AtomSymbol::Element(Element::N), true))),
        b'o' => return Some(Ok((AtomSymbol::Element(Element::O), true))),
        b'p' => return Some(Ok((AtomSymbol::Element(Element::P), true))),
        b's' => return Some(Ok((AtomSymbol::Element(Element::S), true))),
        _ => return None,
    };
    Some(Ok((AtomSymbol::Element(element), false)))
}

#[inline]
fn parse_ascii_element(bytes: &[u8]) -> Option<Element> {
    let symbol = from_utf8(bytes).ok()?;
    Element::from_str(symbol).ok()
}

#[inline]
fn try_element_from_first(
    stream: &mut TokenIter<'_>,
    byte_1: u8,
) -> Result<(AtomSymbol, bool), SmilesError> {
    if byte_1 == b'*' {
        return Ok((AtomSymbol::WildCard, false));
    }
    if !byte_1.is_ascii_alphabetic() {
        return Err(SmilesError::MissingElement);
    }

    let is_aromatic_candidate = byte_1.is_ascii_lowercase();

    if let Some(byte_2) = stream.peek_byte()
        && byte_2.is_ascii_alphabetic()
    {
        if is_aromatic_candidate && byte_2.is_ascii_lowercase() {
            let candidate = [byte_1.to_ascii_uppercase(), byte_2];
            if let Some(element) = parse_ascii_element(&candidate) {
                let _ = stream.next_byte();
                let aromatic = aromatic_from_element(stream.in_bracket, element)?;
                return Ok((AtomSymbol::Element(element), aromatic));
            }
        }
        if !is_aromatic_candidate && byte_2.is_ascii_lowercase() {
            let candidate = [byte_1, byte_2];
            if let Some(element) = parse_ascii_element(&candidate) {
                let _ = stream.next_byte();
                return Ok((AtomSymbol::Element(element), false));
            }
        }
    }

    let one = [if is_aromatic_candidate { byte_1.to_ascii_uppercase() } else { byte_1 }];
    if let Some(element) = parse_ascii_element(&one) {
        let aromatic = if is_aromatic_candidate {
            aromatic_from_element(stream.in_bracket, element)?
        } else {
            false
        };
        return Ok((AtomSymbol::Element(element), aromatic));
    }

    Err(SmilesError::InvalidElementName(char::from(byte_1)))
}

// B, C, N, O, P, S, F, Cl, Br, I,
#[inline]
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
    }
}

#[inline]
fn try_chirality(stream: &mut TokenIter<'_>) -> Result<Option<Chirality>, SmilesError> {
    if stream.peek_byte() != Some(b'@') {
        return Ok(None);
    }
    let _ = stream.next_byte();
    let byte_2 = stream.peek_byte().ok_or(SmilesError::UnexpectedEndOfString)?;
    let chirality = match byte_2 {
        b'@' => {
            let _ = stream.next_byte();
            Chirality::AtAt
        }
        b'T' => {
            let _ = stream.next_byte();
            match stream.peek_byte().ok_or(SmilesError::UnexpectedEndOfString)? {
                b'H' => {
                    let _ = stream.next_byte();
                    let num =
                        try_fold_number::<u8, 1>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_th(num)?
                }
                b'B' => {
                    let _ = stream.next_byte();
                    let num =
                        try_fold_number::<u8, 2>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_tb(num)?
                }
                _ => return Err(SmilesError::InvalidChirality),
            }
        }
        b'A' => {
            let _ = stream.next_byte();
            match stream.peek_byte().ok_or(SmilesError::UnexpectedEndOfString)? {
                b'L' => {
                    let _ = stream.next_byte();
                    let num =
                        try_fold_number::<u8, 1>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_al(num)?
                }
                _ => return Err(SmilesError::InvalidChirality),
            }
        }
        b'S' => {
            let _ = stream.next_byte();
            match stream.peek_byte().ok_or(SmilesError::UnexpectedEndOfString)? {
                b'P' => {
                    let _ = stream.next_byte();
                    let num =
                        try_fold_number::<u8, 1>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_sp(num)?
                }
                _ => return Err(SmilesError::InvalidChirality),
            }
        }
        b'O' => {
            let _ = stream.next_byte();
            match stream.peek_byte().ok_or(SmilesError::UnexpectedEndOfString)? {
                b'H' => {
                    let _ = stream.next_byte();
                    let num =
                        try_fold_number::<u8, 2>(stream).ok_or(SmilesError::InvalidChirality)??;
                    Chirality::try_oh(num)?
                }
                _ => return Err(SmilesError::InvalidChirality),
            }
        }
        b'H' | b'-' | b'+' | b':' | b']' => Chirality::At,
        _ => return Err(SmilesError::InvalidChirality),
    };
    Ok(Some(chirality))
}

#[inline]
fn try_fold_number<B, const MAX_DIGITS: usize>(
    stream: &mut TokenIter<'_>,
) -> Option<Result<B, SmilesError>>
where
    B: TryFrom<u16>,
{
    let mut amount: u16 = 0;
    let mut digits_found = 0;

    let bytes = stream.bytes;
    let len = stream.len;
    let mut position = stream.position;

    while position < len && digits_found < MAX_DIGITS {
        let byte = bytes[position];
        if !byte.is_ascii_digit() {
            break;
        }
        digits_found += 1;
        position += 1;
        let digit = u16::from(byte - b'0');
        match amount.checked_mul(10).and_then(|x| x.checked_add(digit)) {
            Some(val) => amount = val,
            None => return Some(Err(SmilesError::IntegerOverflow)),
        }
    }

    if digits_found == 0 {
        return None;
    }

    stream.position = position;
    Some(B::try_from(amount).map_err(|_| SmilesError::IntegerOverflow))
}

#[inline]
fn hydrogen_count(stream: &mut TokenIter<'_>) -> Result<u8, SmilesError> {
    if stream.peek_byte() == Some(b'H') {
        let _ = stream.next_byte();
        match try_fold_number::<u8, 3>(stream) {
            Some(h) => Ok(h?),
            None => Ok(1),
        }
    } else {
        Ok(0)
    }
}

#[inline]
fn try_charge(stream: &mut TokenIter<'_>) -> Result<Charge, SmilesError> {
    match stream.peek_byte() {
        Some(b'-') => {
            stream.position += 1;
            match stream.peek_byte() {
                Some(b'-') => repeated_charge(stream, b'-'),
                _ => {
                    if let Some(possible_num) = try_fold_number::<i8, 2>(stream) {
                        Charge::try_new(-possible_num?)
                    } else {
                        Charge::try_new(-1)
                    }
                }
            }
        }
        Some(b'+') => {
            stream.position += 1;
            match stream.peek_byte() {
                Some(b'+') => repeated_charge(stream, b'+'),
                _ => {
                    if let Some(possible_num) = try_fold_number::<i8, 2>(stream) {
                        Charge::try_new(possible_num?)
                    } else {
                        Charge::try_new(1)
                    }
                }
            }
        }
        Some(byte) if byte.is_ascii_digit() => {
            if let Some(possible_num) = try_fold_number::<i8, 2>(stream) {
                let magnitude = possible_num?;
                match stream.peek_byte() {
                    Some(b'-') => {
                        stream.position += 1;
                        Charge::try_new(-magnitude)
                    }
                    Some(b'+') => {
                        stream.position += 1;
                        Charge::try_new(magnitude)
                    }
                    _ => Err(SmilesError::UnexpectedCharacter(char::from(byte))),
                }
            } else {
                Err(SmilesError::UnexpectedCharacter(char::from(byte)))
            }
        }
        _ => Ok(Charge::default()),
    }
}

#[inline]
fn repeated_charge(stream: &mut TokenIter<'_>, sign: u8) -> Result<Charge, SmilesError> {
    let mut magnitude: usize = 1;
    while stream.peek_byte() == Some(sign) {
        magnitude = magnitude.checked_add(1).ok_or(if sign == b'+' {
            SmilesError::ChargeOverflow(i8::MAX)
        } else {
            SmilesError::ChargeUnderflow(i8::MIN)
        })?;
        stream.position += 1;
    }

    let magnitude = i8::try_from(magnitude).map_err(|_| {
        if sign == b'+' {
            SmilesError::ChargeOverflow(i8::MAX)
        } else {
            SmilesError::ChargeUnderflow(i8::MIN)
        }
    })?;

    let charge = if sign == b'+' { magnitude } else { -magnitude };
    Charge::try_new(charge)
}

#[inline]
fn try_class(stream: &mut TokenIter<'_>) -> Result<u16, SmilesError> {
    match stream.peek_byte() {
        Some(b':') => {
            stream.position += 1;
            if let Some(possible_num) = try_fold_number::<u16, 3>(stream) {
                possible_num
            } else {
                Err(SmilesError::InvalidClass)
            }
        }
        _ => Ok(0),
    }
}

#[inline]
fn try_bond(byte: u8, bracket: bool) -> Result<Token, SmilesError> {
    let bond = match byte {
        b'-' => {
            if bracket {
                return Err(SmilesError::UnexpectedDash);
            }
            Token::Bond(Bond::Single)
        }
        b'=' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Double));
            }
            Token::Bond(Bond::Double)
        }
        b'#' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Triple));
            }
            Token::Bond(Bond::Triple)
        }
        b'$' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Quadruple));
            }
            Token::Bond(Bond::Quadruple)
        }
        b':' => {
            if bracket {
                return Err(SmilesError::UnexpectedColon);
            }
            Token::Bond(Bond::Aromatic)
        }
        b'/' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Up));
            }
            Token::Bond(Bond::Up)
        }
        b'\\' => {
            if bracket {
                return Err(SmilesError::BondInBracket(Bond::Down));
            }
            Token::Bond(Bond::Down)
        }
        _ => return Err(SmilesError::UnexpectedCharacter(char::from(byte))),
    };
    Ok(bond)
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use elements_rs::Element;

    use super::*;
    use crate::{
        atom::{
            atom_symbol::AtomSymbol,
            bracketed::{charge::Charge, chirality::Chirality},
        },
        bond::{Bond, ring_num::RingNum},
        errors::SmilesError,
        token::Token,
    };

    fn next_ok(input: &str) -> TokenWithSpan {
        TokenIter::from(input).next().expect("expected one token").expect("expected token ok")
    }

    fn next_err(input: &str) -> SmilesErrorWithSpan {
        TokenIter::from(input)
            .next()
            .expect("expected one token")
            .expect_err("expected token error")
    }

    #[test]
    fn parse_token_direct_bracket_state_errors() {
        let mut iter = TokenIter::from(".");
        iter.in_bracket = true;
        assert_eq!(iter.parse_token(b'.'), Err(SmilesError::NonBondInBracket));

        let mut iter = TokenIter::from("[");
        iter.in_bracket = true;
        assert_eq!(iter.parse_token(b'['), Err(SmilesError::UnexpectedLeftBracket));

        let mut iter = TokenIter::from("C");
        iter.in_bracket = true;
        assert_eq!(iter.parse_token(b'C'), Err(SmilesError::UnexpectedBracketedState));

        let mut iter = TokenIter::from("%");
        iter.in_bracket = true;
        assert_eq!(iter.parse_token(b'%'), Err(SmilesError::UnexpectedPercent));

        let mut iter = TokenIter::from("(");
        iter.in_bracket = true;
        assert_eq!(iter.parse_token(b'('), Err(SmilesError::UnexpectedBracketedState));

        let mut iter = TokenIter::from(")");
        iter.in_bracket = true;
        assert_eq!(iter.parse_token(b')'), Err(SmilesError::UnexpectedBracketedState));
    }

    #[test]
    fn parse_token_bracket_atom_covers_isotope_chirality_hydrogens_charge_and_class() {
        let token = next_ok("[13C@H2+2:12]");
        assert_eq!(token.start(), 0);
        assert_eq!(token.end(), "[13C@H2+2:12]".len());
        assert!(matches!(token.token(), Token::Atom(atom) if atom.is_bracket_atom()));
    }

    #[test]
    fn parse_token_unbracketed_atom_errors() {
        let err = next_err("Ac");
        assert_eq!(
            err.smiles_error(),
            SmilesError::InvalidUnbracketedAtom(AtomSymbol::Element(Element::Ac))
        );
        assert_eq!(err.start(), 0);
        assert_eq!(err.end(), 2);
        assert_eq!(err.span().start, 0);
        assert_eq!(err.span().end, 2);

        let err = next_err("q");
        assert_eq!(err.smiles_error(), SmilesError::InvalidElementName('q'));
        assert_eq!(err.start(), 0);
        assert_eq!(err.end(), 1);
        assert_eq!(err.span().start, 0);
        assert_eq!(err.span().end, 1);
    }

    #[test]
    fn parse_token_ring_number_errors() {
        let err = next_err("%");
        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
        assert_eq!(err.start(), 0);
        assert_eq!(err.end(), 1);
        assert_eq!(err.span().start, 0);
        assert_eq!(err.span().end, 1);

        let err = next_err("%9");
        assert_eq!(err.smiles_error(), SmilesError::InvalidRingNumber);
        assert_eq!(err.start(), 0);
        assert_eq!(err.end(), 2);
        assert_eq!(err.span().start, 0);
        assert_eq!(err.span().end, 2);
    }

    #[test]
    fn parse_token_single_digit_ring_closure_success() {
        let token = next_ok("1");
        assert_eq!(token.token(), Token::RingClosure(RingNum::try_new(1).unwrap()));
        assert_eq!(token.span(), 0..1);
    }

    #[test]
    fn aromatic_from_element_branches() {
        assert_eq!(aromatic_from_element(false, Element::C), Ok(true));
        assert_eq!(aromatic_from_element(true, Element::Se), Ok(true));

        assert_eq!(
            aromatic_from_element(false, Element::Se),
            Err(SmilesError::InvalidAromaticElement(Element::Se))
        );

        assert_eq!(
            aromatic_from_element(true, Element::Cl),
            Err(SmilesError::InvalidAromaticElement(Element::Cl))
        );
    }

    #[test]
    fn try_element_from_first_branches() {
        let mut stream = TokenIter::from("");
        assert_eq!(try_element_from_first(&mut stream, b'*'), Ok((AtomSymbol::WildCard, false)));

        let mut stream = TokenIter::from("");
        assert_eq!(
            try_element_from_first(&mut stream, b'c'),
            Ok((AtomSymbol::Element(Element::C), true))
        );

        let mut stream = TokenIter::from("l");
        assert_eq!(
            try_element_from_first(&mut stream, b'C'),
            Ok((AtomSymbol::Element(Element::Cl), false))
        );

        let mut stream = TokenIter::from("e");
        stream.in_bracket = true;
        assert_eq!(
            try_element_from_first(&mut stream, b's'),
            Ok((AtomSymbol::Element(Element::Se), true))
        );

        let mut stream = TokenIter::from("");
        assert_eq!(try_element_from_first(&mut stream, b'1'), Err(SmilesError::MissingElement));

        let mut stream = TokenIter::from("");
        assert_eq!(
            try_element_from_first(&mut stream, b'q'),
            Err(SmilesError::InvalidElementName('q'))
        );
    }

    #[test]
    fn valid_unbracketed_branches() {
        assert!(valid_unbracketed(AtomSymbol::Element(Element::B)));
        assert!(valid_unbracketed(AtomSymbol::Element(Element::C)));
        assert!(valid_unbracketed(AtomSymbol::Element(Element::F)));
        assert!(valid_unbracketed(AtomSymbol::Element(Element::Cl)));
        assert!(valid_unbracketed(AtomSymbol::Element(Element::Br)));
        assert!(valid_unbracketed(AtomSymbol::Element(Element::I)));
        assert!(valid_unbracketed(AtomSymbol::WildCard));
        assert!(!valid_unbracketed(AtomSymbol::Element(Element::Ac)));
    }

    #[test]
    fn try_chirality_branches() {
        let mut stream = TokenIter::from("H");
        assert_eq!(try_chirality(&mut stream), Ok(None));

        let mut stream = TokenIter::from("@@");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::AtAt)));

        let mut stream = TokenIter::from("@]");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::At)));

        let mut stream = TokenIter::from("@H");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::At)));

        let mut stream = TokenIter::from("@+");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::At)));

        let mut stream = TokenIter::from("@AL1");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::try_al(1).unwrap())));

        let mut stream = TokenIter::from("@SP1");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::try_sp(1).unwrap())));

        let mut stream = TokenIter::from("@OH12");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::try_oh(12).unwrap())));

        let mut stream = TokenIter::from("@");
        assert_eq!(try_chirality(&mut stream), Err(SmilesError::UnexpectedEndOfString));

        let mut stream = TokenIter::from("@Z");
        assert_eq!(try_chirality(&mut stream), Err(SmilesError::InvalidChirality));

        let mut stream = TokenIter::from("@AQ");
        assert_eq!(try_chirality(&mut stream), Err(SmilesError::InvalidChirality));

        let mut stream = TokenIter::from("@SQ");
        assert_eq!(try_chirality(&mut stream), Err(SmilesError::InvalidChirality));

        let mut stream = TokenIter::from("@OQ");
        assert_eq!(try_chirality(&mut stream), Err(SmilesError::InvalidChirality));
    }

    #[test]
    fn parse_token_direct_dispatch_covers_remaining_ascii_variants() {
        let mut dot = TokenIter::from(".");
        dot.position = 1;
        assert_eq!(dot.parse_token(b'.').unwrap(), Token::NonBond);

        let mut bracket = TokenIter::from("[C]");
        bracket.position = 1;
        assert_eq!(
            bracket.parse_token(b'[').unwrap(),
            Token::Atom(Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).build())
        );

        let mut invalid = TokenIter::from("Ac");
        invalid.position = 1;
        assert_eq!(
            invalid.parse_token(b'A'),
            Err(SmilesError::InvalidUnbracketedAtom(AtomSymbol::Element(Element::Ac)))
        );

        let mut bond = TokenIter::from("-");
        bond.position = 1;
        assert_eq!(bond.parse_token(b'-').unwrap(), Token::Bond(Bond::Single));

        let mut left = TokenIter::from("(");
        left.position = 1;
        assert_eq!(left.parse_token(b'(').unwrap(), Token::LeftParentheses);

        let mut right = TokenIter::from(")");
        right.position = 1;
        assert_eq!(right.parse_token(b')').unwrap(), Token::RightParentheses);
    }

    #[test]
    fn helper_tables_cover_remaining_allowed_element_paths() {
        for element in [
            Element::B,
            Element::C,
            Element::N,
            Element::O,
            Element::P,
            Element::S,
            Element::Se,
            Element::As,
        ] {
            assert_eq!(aromatic_from_element(true, element), Ok(true));
        }
        for element in [Element::B, Element::C, Element::N, Element::O, Element::S, Element::P] {
            assert_eq!(aromatic_from_element(false, element), Ok(true));
        }
        assert_eq!(
            aromatic_from_element(true, Element::F),
            Err(SmilesError::InvalidAromaticElement(Element::F))
        );
        assert_eq!(
            aromatic_from_element(false, Element::Se),
            Err(SmilesError::InvalidAromaticElement(Element::Se))
        );

        let mut b_only = TokenIter::from("B");
        b_only.position = 1;
        assert_eq!(
            try_organic_subset_from_first(&mut b_only, b'B').unwrap().unwrap(),
            (AtomSymbol::Element(Element::B), false)
        );

        let mut c_only = TokenIter::from("C");
        c_only.position = 1;
        assert_eq!(
            try_organic_subset_from_first(&mut c_only, b'C').unwrap().unwrap(),
            (AtomSymbol::Element(Element::C), false)
        );
    }

    #[test]
    fn try_fold_number_branches_including_overflow() {
        let mut stream = TokenIter::from("123x");
        assert_eq!(try_fold_number::<u16, 3>(&mut stream), Some(Ok(123)));

        let mut stream = TokenIter::from("x");
        assert_eq!(try_fold_number::<u16, 3>(&mut stream), None);

        let mut stream = TokenIter::from("700000");
        assert_eq!(try_fold_number::<u16, 6>(&mut stream), Some(Err(SmilesError::IntegerOverflow)));

        let mut stream = TokenIter::from("300");
        assert_eq!(try_fold_number::<u8, 3>(&mut stream), Some(Err(SmilesError::IntegerOverflow)));
    }

    #[test]
    fn hydrogen_count_branches() {
        let mut stream = TokenIter::from("H2");
        assert_eq!(hydrogen_count(&mut stream), Ok(2));

        let mut stream = TokenIter::from("H");
        assert_eq!(hydrogen_count(&mut stream), Ok(1));

        let mut stream = TokenIter::from("C");
        assert_eq!(hydrogen_count(&mut stream), Ok(0));
    }

    #[test]
    fn try_charge_branches() {
        let mut stream = TokenIter::from("-");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(-1).unwrap()));

        let mut stream = TokenIter::from("--");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(-2).unwrap()));

        let mut stream = TokenIter::from("-15");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(-15).unwrap()));

        let mut stream = TokenIter::from("-16");
        assert_eq!(try_charge(&mut stream), Err(SmilesError::ChargeUnderflow(-16)));

        let mut stream = TokenIter::from("+");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(1).unwrap()));

        let mut stream = TokenIter::from("++");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(2).unwrap()));

        let mut stream = TokenIter::from("+15");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(15).unwrap()));

        let mut stream = TokenIter::from("+16");
        assert_eq!(try_charge(&mut stream), Err(SmilesError::ChargeOverflow(16)));

        let mut stream = TokenIter::from("C");
        assert_eq!(try_charge(&mut stream), Ok(Charge::default()));
    }

    #[test]
    fn try_class_branches() {
        let mut stream = TokenIter::from(":12");
        assert_eq!(try_class(&mut stream), Ok(12));

        let mut stream = TokenIter::from(":");
        assert_eq!(try_class(&mut stream), Err(SmilesError::InvalidClass));

        let mut stream = TokenIter::from("C");
        assert_eq!(try_class(&mut stream), Ok(0));
    }

    #[test]
    fn try_bond_branches() {
        let ok_cases = [
            ('-', Token::Bond(Bond::Single)),
            ('=', Token::Bond(Bond::Double)),
            ('#', Token::Bond(Bond::Triple)),
            ('$', Token::Bond(Bond::Quadruple)),
            (':', Token::Bond(Bond::Aromatic)),
            ('/', Token::Bond(Bond::Up)),
            ('\\', Token::Bond(Bond::Down)),
        ];

        for (ch, expected) in ok_cases {
            assert_eq!(try_bond(ch as u8, false), Ok(expected));
        }

        assert_eq!(try_bond(b'-', true), Err(SmilesError::UnexpectedDash));
        assert_eq!(try_bond(b':', true), Err(SmilesError::UnexpectedColon));
        assert_eq!(try_bond(b'=', true), Err(SmilesError::BondInBracket(Bond::Double)));
        assert_eq!(try_bond(b'#', true), Err(SmilesError::BondInBracket(Bond::Triple)));
        assert_eq!(try_bond(b'$', true), Err(SmilesError::BondInBracket(Bond::Quadruple)));
        assert_eq!(try_bond(b'/', true), Err(SmilesError::BondInBracket(Bond::Up)));
        assert_eq!(try_bond(b'\\', true), Err(SmilesError::BondInBracket(Bond::Down)));
        assert_eq!(try_bond(b'x', false), Err(SmilesError::UnexpectedCharacter('x')));
    }

    #[test]
    fn iterator_error_span_mapping_smoke_test() {
        let err = next_err("Ac");
        assert_eq!(
            err.smiles_error(),
            SmilesError::InvalidUnbracketedAtom(AtomSymbol::Element(Element::Ac))
        );
        assert_eq!(err.start(), 0);
        assert_eq!(err.end(), 2);
    }

    #[test]
    fn try_chirality_th_form_should_parse() {
        let mut stream = TokenIter::from("@TH1");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::try_th(1).unwrap())));
    }

    #[test]
    fn try_chirality_tb_form_should_parse() {
        let mut stream = TokenIter::from("@TB10");
        assert_eq!(try_chirality(&mut stream), Ok(Some(Chirality::try_tb(10).unwrap())));
    }
    #[test]
    fn token_iter_parses_bracket_atom_with_th_chirality() {
        let token = next_ok("[C@TH1]");
        assert!(matches!(token.token(), Token::Atom(atom) if atom.is_bracket_atom()));
    }
    #[test]
    fn test_charge_parsing() {
        // Define test cases with expected charge parsing behavior
        let test_cases = vec![
            ("[Fe+++]", 3), // Fe with +3 charge
            ("[Fe-1]", -1), // Fe with -1 charge
            ("[Fe+]", 1),   // Fe with +1 charge
            ("[Fe2+]", 2),  // Fe with +2 charge
            ("[Fe-2]", -2), // Fe with -2 charge
        ];

        for (smiles, expected_charge) in test_cases {
            let tokens: Vec<Result<TokenWithSpan, SmilesErrorWithSpan>> =
                TokenIter::from(smiles).collect();
            for token in tokens {
                match token {
                    Ok(token_with_span) => {
                        let possible_atom = token_with_span.token();
                        match possible_atom {
                            Token::Atom(atom) => {
                                assert_eq!(atom.charge_value(), expected_charge);
                            }
                            _ => panic!("Token: {possible_atom:?} should be an atom token!"),
                        }
                    }
                    Err(e) => panic!("{e} in {smiles}"),
                }
            }
        }
    }
    #[test]
    fn try_charge_supports_number_before_sign() {
        let mut stream = TokenIter::from("2-");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(-2).unwrap()));

        let mut stream = TokenIter::from("2+");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(2).unwrap()));

        let mut stream = TokenIter::from("15-");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(-15).unwrap()));

        let mut stream = TokenIter::from("15+");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(15).unwrap()));
    }
    #[test]
    fn test_charge_parsing_number_before_sign() {
        let test_cases = vec![("[C2-]", -2), ("[C2+]", 2), ("[Fe15+]", 15), ("[Fe15-]", -15)];

        for (smiles, expected_charge) in test_cases {
            let tokens: Vec<_> = TokenIter::from(smiles).collect();
            for token in tokens {
                match token {
                    Ok(token_with_span) => {
                        match token_with_span.token() {
                            Token::Atom(atom) => {
                                assert_eq!(atom.charge_value(), expected_charge);
                            }
                            other => panic!("Token {other:?} should be an atom token!"),
                        }
                    }
                    Err(e) => panic!("{e} in {smiles}"),
                }
            }
        }
    }

    #[test]
    fn parse_token_misc_ascii_and_unicode_error_paths() {
        let token = next_ok(".");
        assert_eq!(token.token(), Token::NonBond);

        let token = next_ok("(");
        assert_eq!(token.token(), Token::LeftParentheses);

        let token = next_ok(")");
        assert_eq!(token.token(), Token::RightParentheses);

        let err = next_err("€");
        assert_eq!(err.smiles_error(), SmilesError::UnexpectedUnicodeCharacter);
        assert_eq!(err.start(), 0);
        assert_eq!(err.end(), 3);
    }

    #[test]
    fn utf8_char_width_covers_all_width_buckets() {
        assert_eq!(utf8_char_width(b'a'), 1);
        assert_eq!(utf8_char_width(0xc2), 2);
        assert_eq!(utf8_char_width(0xe2), 3);
        assert_eq!(utf8_char_width(0xf0), 4);
    }

    #[test]
    fn organic_subset_fast_path_covers_remaining_variants() {
        let cases = [
            ("B", AtomSymbol::Element(Element::B), false),
            ("F", AtomSymbol::Element(Element::F), false),
            ("I", AtomSymbol::Element(Element::I), false),
            ("b", AtomSymbol::Element(Element::B), true),
            ("o", AtomSymbol::Element(Element::O), true),
            ("p", AtomSymbol::Element(Element::P), true),
            ("s", AtomSymbol::Element(Element::S), true),
            ("*", AtomSymbol::WildCard, false),
        ];

        for (input, symbol, aromatic) in cases {
            let token = next_ok(input);
            assert_eq!(token.token(), Token::Atom(Atom::new_organic_subset(symbol, aromatic)));
        }
    }

    #[test]
    fn try_element_and_general_element_path_cover_bracket_specific_cases() {
        let mut stream = TokenIter::from("Xe");
        assert_eq!(try_element(&mut stream), Ok((AtomSymbol::Element(Element::Xe), false)));

        let mut stream = TokenIter::from("se");
        stream.in_bracket = true;
        assert_eq!(try_element(&mut stream), Ok((AtomSymbol::Element(Element::Se), true)));

        let mut stream = TokenIter::from("q");
        assert_eq!(try_element(&mut stream), Err(SmilesError::InvalidElementName('q')));

        let mut stream = TokenIter::from("");
        assert_eq!(try_element(&mut stream), Err(SmilesError::MissingElement));
    }

    #[test]
    fn parse_ascii_element_and_valid_unbracketed_cover_false_cases() {
        assert_eq!(parse_ascii_element(b"Xe"), Some(Element::Xe));
        assert_eq!(parse_ascii_element(b"??"), None);
        assert!(!valid_unbracketed(AtomSymbol::Element(Element::Xe)));
    }

    #[test]
    fn try_charge_handles_digit_prefixed_forms_and_errors() {
        let mut stream = TokenIter::from("2+");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(2).unwrap()));

        let mut stream = TokenIter::from("2-");
        assert_eq!(try_charge(&mut stream), Ok(Charge::try_new(-2).unwrap()));

        let mut stream = TokenIter::from("2x");
        assert_eq!(try_charge(&mut stream), Err(SmilesError::UnexpectedCharacter('2')));
    }

    #[test]
    fn try_class_parses_large_values_and_try_bond_covers_remaining_paths() {
        let mut stream = TokenIter::from(":999");
        assert_eq!(try_class(&mut stream), Ok(999));

        assert_eq!(try_bond(b'=', false), Ok(Token::Bond(Bond::Double)));
        assert_eq!(try_bond(b'#', false), Ok(Token::Bond(Bond::Triple)));
        assert_eq!(try_bond(b'$', false), Ok(Token::Bond(Bond::Quadruple)));
        assert_eq!(try_bond(b'/', false), Ok(Token::Bond(Bond::Up)));
        assert_eq!(try_bond(b'\\', false), Ok(Token::Bond(Bond::Down)));
    }

    #[test]
    fn repeated_positive_charge_overflow_should_error_instead_of_panicking() {
        let smiles = format!("[n{}]", "+".repeat(200));

        let result = std::panic::catch_unwind(|| TokenIter::from(smiles.as_str()).next());
        assert!(result.is_ok(), "tokenization panicked on oversized repeated positive charge");

        let token = result.unwrap().expect("expected one token");
        assert!(
            matches!(
                token,
                Err(ref err) if matches!(err.smiles_error(), SmilesError::ChargeOverflow(_))
            ),
            "expected ChargeOverflow error, got {token:?}"
        );
    }

    #[test]
    fn repeated_negative_charge_underflow_should_error_instead_of_panicking() {
        let smiles = format!("[n{}]", "-".repeat(200));

        let result = std::panic::catch_unwind(|| TokenIter::from(smiles.as_str()).next());
        assert!(result.is_ok(), "tokenization panicked on oversized repeated negative charge");

        let token = result.unwrap().expect("expected one token");
        assert!(
            matches!(
                token,
                Err(ref err) if matches!(err.smiles_error(), SmilesError::ChargeUnderflow(_))
            ),
            "expected ChargeUnderflow error, got {token:?}"
        );
    }
}
