//! Defines errors used in the SMILES parser.

use std::num::TryFromIntError;

use elements_rs::Element;

use crate::token::Bond;

/// The errors that could occur during SMILES parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesError {
    /// Missing Element
    MissingElement,
    /// Invalid Isotope value passed
    InvalidIsotope,
    /// Error indicating invalid Element name
    InvalidElementName(char),
    /// Error indicating that an invalid number was encountered.
    InvalidNumber,
    /// Error indicating that an unexpected character was encountered.
    UnexpectedCharacter {
        /// The unexpected character.
        character: char,
    },
    /// An unexpected left bracket `[` was found
    UnexpectedLeftBracket,
    /// An unexpected right bracket `]` was found
    UnexpectedRightBracket,
    /// A closing `]` bracket was not found
    UnclosedBracket,
    /// A non bare element found outside of brackets
    ElementRequiresBrackets,
    /// found `[..]` that did not contain an element
    MissingBracketElement,
    /// Element forbidden to be written as aromatic here
    InvalidAromaticElement {
        /// The forbidden aromatic element
        element: Element,
    },
    /// Integer Overflow
    IntegerOverflow,
    /// Ring Number Overflow (greater than 99)
    RingNumberOverflow(u8),
    /// A charge is below allowed minimum (-15)
    ChargeUnderflow(i8),
    /// A charge is over the allowed maximum (15)
    ChargeOverflow(i8),
    /// Bond Inside Bracket
    BondInBracket(Bond),
    /// Non Bond in Bracket
    NonBondInBracket,
    /// Wrapper for `element_rs` errors
    ElementsRs(elements_rs::errors::Error),
}

impl From<elements_rs::errors::Error> for SmilesError {
    fn from(e: elements_rs::errors::Error) -> Self {
        SmilesError::ElementsRs(e)
    }
}

impl From<TryFromIntError> for SmilesError {
    fn from(_: TryFromIntError) -> Self {
        SmilesError::InvalidNumber
    }
}

