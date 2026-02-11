//! Submodule creating the `TokenIter` struct, which is an iterator over
//! the `Token`s found in a provided string.

use std::str::FromStr;

use elements_rs::Element;

use crate::{
    errors::SmilesError,
    token::{AtomSymbol, BracketedAtom, Token, UnbracketedAtom},
};

/// An iterator over the tokens found in a SMILES string.
pub struct TokenIter<'a> {
    /// The peekable chars iterator
    chars: std::iter::Peekable<std::str::Chars<'a>>,
    /// Denotes whether currently inside brackets
    in_bracket: bool,
    /// the previous char
    prev_char: Option<char>,
}

impl<'a> From<&'a str> for TokenIter<'a> {
    fn from(s: &'a str) -> Self {
        TokenIter { chars: s.chars().peekable(), in_bracket: false, prev_char: None }
    }
}

impl TokenIter<'_> {
    fn parse_token(&mut self, current_char: char) -> Result<Token, SmilesError> {
        let token = match current_char {
            '.' => {
                if !self.in_bracket {
                    Token::NonBond
                } else {
                    return Err(SmilesError::NonBondInBracket);
                }
            }
            '[' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedLeftBracket);
                } else {
                    self.in_bracket = true;
                    let possible_bracket_atom = BracketedAtom::builder();
                    let try_isotope = try_fold_number(self);
                    if let Some(isotope) = try_isotope {
                        let val = isotope?;
                        possible_bracket_atom.with_isotope(val);
                    }

                    let bracket_atom = possible_bracket_atom.build();
                    Token::BracketedAtom(bracket_atom)
                }
            }
            _ => return Err(SmilesError::UnexpectedCharacter { character: current_char }),
        };
        Ok(token)
    }
}

impl Iterator for TokenIter<'_> {
    type Item = Result<Token, SmilesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.chars.next().map(|current_char| self.parse_token(current_char))
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
    if allowed { Ok(true) } else { Err(SmilesError::InvalidAromaticElement { element }) }
}

fn try_element(stream: &mut TokenIter<'_>) -> Result<AtomSymbol, SmilesError> {
    if matches!(stream.chars.peek(), Some('*')) {
        stream.chars.next();
        return Ok(AtomSymbol::default());
    }
    let Some(&char_1) = stream.chars.peek() else {
        return Err(SmilesError::MissingElement);
    };
    if !char_1.is_alphabetic() {
        return Err(SmilesError::MissingElement);
    }
    let try_candidate = |val: &str| -> Option<Element> { Element::from_str(val).ok() };
    stream.chars.next();
    let first_char = char_1.to_string();
    if let Some(&char_2) = stream.chars.peek() {
        if char_2.is_alphabetic() {
            let two_chars = format!("{char_1}{char_2}");
            if let Some(element) = try_candidate(&two_chars) {
                stream.chars.next();
                return Ok(AtomSymbol::Element(element));
            }
        }
    }
    if let Some(element) = try_candidate(&first_char) {
        return Ok(AtomSymbol::Element(element));
    }
    Err(SmilesError::InvalidElementName(char_1))
}

fn try_fold_number<B>(stream: &mut TokenIter<'_>) -> Option<Result<B, SmilesError>>
where
    B: TryFrom<u32>,
{
    let mut seen_any = false;
    let mut amount: u32 = 0;

    while let Some(char) = stream.chars.peek() {
        let digit = match char.to_digit(10) {
            Some(d) => d,
            None => break,
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
