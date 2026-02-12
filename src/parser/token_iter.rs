//! Submodule creating the `TokenIter` struct, which is an iterator over
//! the `Token`s found in a provided string.

use std::str::FromStr;

use elements_rs::Element;

use crate::{
    errors::SmilesError,
    token::{AtomSymbol, BracketedAtom, Chirality, Token, UnbracketedAtom},
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
                    if let Some(isotope) = try_fold_number(self) {
                        possible_bracket_atom.with_isotope(isotope?);
                    }
                    let (atom, aromatic) = try_element(self)?;
                    possible_bracket_atom.with_element(atom.element()).with_aromatic(aromatic);
                    todo!("add chirality method, hydrogen count, charge, and class");
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

fn try_element(stream: &mut TokenIter<'_>) -> Result<(AtomSymbol, bool), SmilesError> {
    if matches!(stream.chars.peek(), Some('*')) {
        stream.chars.next();
        return Ok((AtomSymbol::WildCard, false));
    }
    let Some(&char_1) = stream.chars.peek() else {
        return Err(SmilesError::MissingElement);
    };
    if !char_1.is_alphabetic() {
        return Err(SmilesError::MissingElement);
    }
    let is_aromatic_candidate = char_1.is_ascii_lowercase();
    stream.chars.next();

    let try_candidate = |val: &str| -> Option<Element> { Element::from_str(val).ok() };

    if let Some(&char_2) = stream.chars.peek() {
        if char_2.is_alphabetic() {
            if is_aromatic_candidate && char_2.is_ascii_lowercase() {
                let candidate = format!("{}{}", char_1.to_ascii_uppercase(), char_2);
                if let Some(element) = try_candidate(&candidate) {
                    stream.chars.next();
                    let aromatic = aromatic_from_element(stream.in_bracket, element)?;
                    return Ok((AtomSymbol::Element(element), aromatic));
                }
            }
            if !is_aromatic_candidate && char_2.is_ascii_lowercase() {
                let candidate = format!("{}{}", char_1, char_2);
                if let Some(element) = try_candidate(&candidate) {
                    stream.chars.next();
                    return Ok((AtomSymbol::Element(element), false));
                }
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

fn try_chirality(stream: &mut TokenIter<'_>) -> Result<Option<Chirality>,SmilesError> {
    let chirality = if let Some(char_1) = stream.chars.peek() {
        stream.chars.next();
        if let Some(char_2) = stream.chars.peek() {
            match (char_1, char_2) {
                ('@','@') => {
                    stream.chars.next();
                    Ok(Some(Chirality::AtAt))
                },
                ('@', 'T') => {
                    stream.chars.next();
                    if let Some(char_3) = stream.chars.peek() {
                        match char_3 {
                            'H' => {
                                stream.chars.next();
                                let possible_num = try_fold_number(stream);
                                match possible_num {
                                    None => Err(SmilesError::InvalidChirality),
                                    Some(num) => Chirality::try_th(num?)
                                }
                            },
                            'B' => {
                                stream.chars.next();
                                let possible_num = try_fold_number(stream);
                                match possible_num {
                                    None => Err(SmilesError::InvalidChirality),
                                    Some(num) => Chirality::try_tb(num?)
                                }
                            },
                            _ => Err(SmilesError::InvalidChirality)
                        }
                    } else {
                        Err(SmilesError::InvalidChirality)
                    }
                }
            }
        } else {
            Err(SmilesError::UnexpectedEndOfString)            
        }
    } else {
        Err(SmilesError::UnexpectedEndOfString)
    };
    chirality
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
