//! Submodule creating the `TokenIter` struct, which is an iterator over
//! the `Token`s found in a provided string.

use elements_rs::
    Element;

use crate::{errors::SmilesError, token::Token};

/// An iterator over the tokens found in a SMILES string.
pub struct TokenIter<'a> {
    /// The peekable chars iterator
    chars: std::iter::Peekable<std::str::Chars<'a>>,
}

impl<'a> From<&'a str> for TokenIter<'a> {
    fn from(s: &'a str) -> Self {
        TokenIter { chars: s.chars().peekable() }
    }
}

impl TokenIter<'_> {
    fn parse_token(&mut self, current_char: char) -> Result<Token, crate::errors::SmilesError> {
        Ok(match current_char {
            '(' => Token::OpenRoundBracket,
            ')' => Token::CloseRoundBracket,
            '=' => Token::Equal,
            '#' => Token::Hashtag,
            '$' => Token::Dollar,
            '.' => Token::Dot,
            ':' => Token::Colon,
            '/' => Token::ForwardSlash,
            '\\' => Token::BackSlash,
            '[' => {
                let molecule_iter = self.chars.by_ref().take_while(|c| *c != ']');
                todo!()
            }
            ']' => Token::RightSquareBracket,
            number @ '0'..='9' => {
                if let Some(num) = number.to_digit(10) {
                    let label = u8::try_from(num)?;
                    return Ok(Token::Label(label));
                } else {
                    return Err(SmilesError::InvalidNumber);
                }
            }
            maybe_element_char => {
                if let Some(next_char) = self.chars.peek()
                    && let Ok(element) = Element::try_from([maybe_element_char, *next_char])
                {
                    self.chars.next();
                    todo!()
                }
                if let Ok(element) = Element::try_from(maybe_element_char) {
                    todo!()
                }
                return Err(SmilesError::UnexpectedCharacter { character: maybe_element_char });
            }
        })
    }
}

impl Iterator for TokenIter<'_> {
    type Item = Result<crate::token::Token, SmilesError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.chars.next().map(|current_char| self.parse_token(current_char))
    }
}
