//! Defines errors used in the SMILES parser.

use core::fmt;
use std::{num::TryFromIntError, ops::Range};

use elements_rs::Element;

use crate::{atom_symbol::AtomSymbol, bond::Bond};

/// The errors that could occur during SMILES parsing.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
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
    UnexpectedCharacter(char),
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
    InvalidAromaticElement(Element),
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
    /// Specified Chirality is not a valid form
    InvalidChirality,
    /// Unexpected end of string
    UnexpectedEndOfString,
    /// The class is not valid
    InvalidClass,
    /// Non organic element found out of bracket
    InvalidUnbracketedAtom(AtomSymbol),
    /// Unexpectedly inside of brackets
    UnexpectedBracketedState,
    /// An unexpected `-` has been found
    UnexpectedDash,
    /// An unexpected `:` has been found
    UnexpectedColon,
    /// An unexpected `%` has been found
    UnexpectedPercent,
    /// An invalid ring number has been found
    InvalidRingNumber,
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
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

/// Wraps the Smiles error adding the location of where the error was found
pub struct SmilesErrorWithSpan {
    /// The [`SmilesError`]
    smiles_error: SmilesError,
    /// The span as `usize`
    span: Range<usize>,
}

impl SmilesErrorWithSpan {
    /// Creates a new error from the [`SmilesError`] and the `span`
    pub fn new(smiles_error: SmilesError, start: usize, end: usize) -> Self {
        Self { smiles_error, span: Range { start: start, end: end } }
    }
    /// Returns the [`SmilesError`]
    pub fn smiles_error(&self) -> SmilesError {
        self.smiles_error
    }
    /// Returns the start of the span
    pub fn start(&self) -> usize {
        self.span.start
    }
    /// Returns the end of the span
    pub fn end(&self) -> usize {
        self.span.end
    }
}

impl fmt::Display for SmilesErrorWithSpan {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.smiles_error() {
            SmilesError::InvalidRingNumber
            | SmilesError::UnexpectedPercent
            | SmilesError::UnexpectedColon
            | SmilesError::UnexpectedDash
            | SmilesError::UnexpectedBracketedState
            | SmilesError::InvalidClass
            | SmilesError::UnexpectedEndOfString
            | SmilesError::InvalidChirality
            | SmilesError::NonBondInBracket
            | SmilesError::IntegerOverflow
            | SmilesError::MissingBracketElement
            | SmilesError::ElementRequiresBrackets
            | SmilesError::UnclosedBracket
            | SmilesError::UnexpectedRightBracket
            | SmilesError::UnexpectedLeftBracket
            | SmilesError::InvalidNumber
            | SmilesError::InvalidIsotope
            | SmilesError::MissingElement => {
                write!(f, "{} at: {} end: {}", self.smiles_error(), self.start(), self.end())
            }

            SmilesError::UnexpectedCharacter(e) | SmilesError::InvalidElementName(e) => {
                write!(f, "{}:{} at: {} end: {}", self.smiles_error(), e, self.start(), self.end())
            }
            SmilesError::InvalidAromaticElement(e) => {
                write!(f, "{}:{} at: {} end: {}", self.smiles_error(), e, self.start(), self.end())
            }

            SmilesError::RingNumberOverflow(r) => {
                write!(f, "{}:{} at: {} end: {}", self.smiles_error(), r, self.start(), self.end())
            }
            SmilesError::ChargeOverflow(c) | SmilesError::ChargeUnderflow(c) => {
                write!(f, "{}:{} at: {} end: {}", self.smiles_error(), c, self.start(), self.end())
            }
            SmilesError::BondInBracket(bond) => {
                write!(
                    f,
                    "{}:{} at: {} end: {}",
                    self.smiles_error(),
                    bond,
                    self.start(),
                    self.end()
                )
            }

            SmilesError::ElementsRs(error) => {
                write!(
                    f,
                    "{}:{} at: {} end: {}",
                    self.smiles_error(),
                    error,
                    self.start(),
                    self.end()
                )
            }

            SmilesError::InvalidUnbracketedAtom(atom_symbol) => {
                write!(
                    f,
                    "{}:{} at: {} end: {}",
                    self.smiles_error(),
                    atom_symbol,
                    self.start(),
                    self.end()
                )
            }
        }
    }
}
