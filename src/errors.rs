//! Defines errors used in the SMILES parser.

use std::num::TryFromIntError;

use elements_rs::Element;

/// The errors that could occur during SMILES parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesError {
    /// Error indicating incomplete Element name
    IncompleteElementName(char),
    /// Error indicating that an unknown element was encountered.
    Element(elements_rs::errors::Error),
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
}

impl From<TryFromIntError> for SmilesError {
    fn from(_: TryFromIntError) -> Self {
        SmilesError::InvalidNumber
    }
}

impl From<elements_rs::errors::Error> for SmilesError {
    fn from(value: elements_rs::errors::Error) -> Self {
        SmilesError::Element(value)
    }
}
