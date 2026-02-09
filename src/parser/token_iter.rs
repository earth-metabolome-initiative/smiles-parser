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
  
        Ok(token)
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

fn try_fold_number<A, B>(stream: &mut TokenIter<'_>) -> Option<Result<B, SmilesError>>
where
    A: Iterator<Item = char>,
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
