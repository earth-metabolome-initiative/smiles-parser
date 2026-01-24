//! Submodule creating the `TokenIter` struct, which is an iterator over
//! the `Token`s found in a provided string.

use elements_rs::Element;

use crate::{errors::SmilesError, token::Token};

/// An iterator over the tokens found in a SMILES string.
pub struct TokenIter<'a> {
    /// The peekable chars iterator
    chars: std::iter::Peekable<std::str::Chars<'a>>,
    /// Denotes whether currently inside brackets
    in_bracket: bool,
}

impl<'a> From<&'a str> for TokenIter<'a> {
    fn from(s: &'a str) -> Self {
        TokenIter { chars: s.chars().peekable(), in_bracket: false }
    }
}

impl TokenIter<'_> {
    fn parse_token(&mut self, current_char: char) -> Result<Token, SmilesError> {
        Ok(match current_char {
            '&' => Token::Ampersand,
            atom_char if atom_char.is_ascii_alphabetic() => {
                if atom_char.is_ascii_uppercase() {
                    let two = self
                        .chars
                        .peek()
                        .copied()
                        .filter(char::is_ascii_lowercase)
                        .and_then(|c| Element::try_from([atom_char, c]).ok());

                    if let Some(element) = two {
                        self.chars.next();
                        Token::Atom { element, aromatic: false }
                    } else if let Ok(element) = Element::try_from(atom_char) {
                        Token::Atom { element, aromatic: false }
                    } else {
                        return Err(SmilesError::UnexpectedCharacter { character: atom_char });
                    }
                } else {
                    let upper = atom_char.to_ascii_uppercase();
                    // There's only two possible two char aromatics
                    let two = match atom_char {
                        'a' => {
                            self.chars
                                .peek()
                                .copied()
                                .filter(|&c| c == 's')
                                .and_then(|c| Element::try_from([upper, c]).ok())
                        }
                        's' => {
                            self.chars
                                .peek()
                                .copied()
                                .filter(|&c| c == 'e')
                                .and_then(|c| Element::try_from([upper, c]).ok())
                        }
                        _ => None,
                    };
                    if let Some(element) = two {
                        self.chars.next();
                        Token::Atom { element, aromatic: self.aromatic_from_element(element)? }
                    } else if let Ok(element) = Element::try_from(upper) {
                        Token::Atom { element, aromatic: self.aromatic_from_element(element)? }
                    } else {
                        return Err(SmilesError::UnexpectedCharacter { character: atom_char });
                    }
                }
            }
            '@' => Token::AtSign,
            '\\' => Token::BackSlash,
            ':' => Token::Colon,
            '$' => Token::Dollar,
            '.' => Token::Dot,
            '=' => Token::Equal,
            '/' => Token::ForwardSlash,
            '#' => Token::Hashtag,
            label if label.is_ascii_digit() => {
                // safe to get value by casting to ascii value and offset.
                Token::Label(label as u8 - b'0')
            }
            '(' => Token::LeftParentheses,
            '[' => {
                if self.in_bracket {
                    // can't have nested brackets no need to parse, already invalid
                    return Err(SmilesError::UnexpectedLeftBracket);
                }
                self.in_bracket = true;
                Token::LeftSquareBracket
            }
            '-' => Token::Minus,
            '%' => Token::Percent,
            '+' => Token::Plus,
            ')' => Token::RightParentheses,
            ']' => {
                if !self.in_bracket {
                    return Err(SmilesError::UnexpectedRightBracket);
                }
                self.in_bracket = false;
                Token::RightSquareBracket
            }
            c => return Err(SmilesError::UnexpectedCharacter { character: c }),
        })
    }

    /// determines whether an aromatic is valid for a given [`Token::Atom`]
    ///
    /// # Parameters
    /// - the tokens as `self`
    /// - the [`Element`] being passed
    fn aromatic_from_element(&self, element: Element) -> Result<bool, SmilesError> {
        let allowed = if self.in_bracket {
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
}

impl Iterator for TokenIter<'_> {
    type Item = Result<Token, SmilesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.chars.next().map(|current_char| self.parse_token(current_char))
    }
}
