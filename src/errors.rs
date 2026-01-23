//! Defines errors used in the SMILES parser.

use std::num::TryFromIntError;

/// The errors that could occur during SMILES parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesError {
    /// Error indicating that an unknown element was encountered.
    Element(elements_rs::errors::Error),
    /// Error indicating that an invalid number was encountered.
    InvalidNumber,
    /// Error indicating that an unexpected character was encountered.
    UnexpectedCharacter {
        /// The unexpected character.
        character: char,
    },
}

impl From<TryFromIntError> for SmilesError {
    fn from(_: TryFromIntError) -> Self {
        SmilesError::InvalidNumber
    }
}
