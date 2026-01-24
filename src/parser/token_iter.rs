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
    fn parse_token(&mut self, current_char: char) -> Result<Token, crate::errors::SmilesError> {
        Ok(match current_char {
            '(' => Token::LeftParentheses,
            ')' => Token::RightParentheses,
            '=' => Token::Equal,
            '#' => Token::Hashtag,
            '$' => Token::Dollar,
            '.' => Token::Dot,
            ':' => Token::Colon,
            '/' => Token::ForwardSlash,
            '\\' => Token::BackSlash,
            '[' => {
                if self.in_bracket {
                    return Err(SmilesError::UnexpectedLeftBracket);
                }
                self.in_bracket = true;
                Token::LeftSquareBracket
            }
            ']' => {
                if !self.in_bracket {
                    return Err(SmilesError::UnexpectedRightBracket);
                }
                self.in_bracket = false;
                Token::RightSquareBracket
            }
            num if num.is_ascii_digit() => {
                // safe to get value by casting to ascii value and offset. 
                Token::Label(num as u8 - b'0')
            }

            element_char if element_char.is_ascii_alphabetic() => {
                if element_char.is_uppercase() {
                    let next = self.chars.peek();
                    match element_char {
                        'A' => {
                            match next {
                                Some(&c) => {
                                    match c {
                                        'c' |'g'|'l'|'m'|'r'|'s'|'t'|'u' => {
                                            let element = Element::try_from([element_char, c])?;
                                            self.chars.next();
                                            Token::Atom { element: element, aromatic: false }
                                        }
                                    }
                                },
                                None => return Err(SmilesError::IncompleteElementName(element_char))
                            }
                        },
                        'B' => {},
                        'C' => {},
                        'D' => {},
                        'E' => {},
                        'F' => {},
                        'G' => {},
                        'H' => {},
                        'I' => {},
                        'J' => {},
                        'K' => {},
                        'L' => {},
                        'M' => {},
                        'N' => {},
                        'O' => {},
                        'P' => {},
                        'Q' => {},
                        'R' => {},
                        'S' => {},
                        'T' => {},
                        'U' => {},
                        'V' => {},
                        'W' => {},
                        'X' => {},
                        'Y' => {},
                        'Z' => {},
                    }
                } else {

                }
            }
            // maybe_element_char => {
            //     let molecule_iter = self.chars.by_ref().take_while(|c| c.is_alphabetic());
            //     if let Some(next_char) = self.chars.peek()
            //         && let Ok(element) = Element::try_from([maybe_element_char, *next_char])
            //     {
            //         self.chars.next();
            //         todo!()
            //     }
            //     if let Ok(element) = Element::try_from(maybe_element_char) {
            //         todo!()
            //     }
            //     return Err(SmilesError::UnexpectedCharacter { character: maybe_element_char });
            // }
            c => return Err(SmilesError::UnexpectedCharacter { character: c }),
        })
    }

    /// determines whether an aromatic is valid for a given [`Token::Atom`]
    /// 
    /// # Parameters
    /// - the tokens as `self`
    /// - the [`Element`] being passed
    fn aromatic_from_element(
        &self,
        element: Element
    ) -> Result<bool, SmilesError> {
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
    type Item = Result<crate::token::Token, SmilesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.chars.next().map(|current_char| self.parse_token(current_char))
    }
}
