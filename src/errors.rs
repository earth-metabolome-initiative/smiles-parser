//! Defines errors used in the SMILES parser.

use core::fmt;
use std::{num::TryFromIntError, ops::Range};

use elements_rs::Element;

use crate::{atom_symbol::AtomSymbol, bond::Bond};

/// The errors that could occur during SMILES parsing.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SmilesError {
    /// Bond Inside Bracket
    BondInBracket(Bond),
    /// A charge is over the allowed maximum (15)
    ChargeOverflow(i8),
    /// A charge is below allowed minimum (-15)
    ChargeUnderflow(i8),
    /// A non bare element found outside of brackets
    ElementRequiresBrackets,
    /// Wrapper for `element_rs` errors
    ElementsRs(elements_rs::errors::Error),
    /// Element forbidden to be written as aromatic here
    InvalidAromaticElement(Element),
    /// Specified Chirality is not a valid form
    InvalidChirality,
    /// The class is not valid
    InvalidClass,
    /// Error indicating invalid Element name
    InvalidElementName(char),
    /// Invalid Isotope value passed
    InvalidIsotope,
    /// Error indicating that an invalid number was encountered.
    InvalidNumber,
    /// Integer Overflow
    IntegerOverflow,
    /// Non organic element found out of bracket
    InvalidUnbracketedAtom(AtomSymbol),
    /// An invalid ring number has been found
    InvalidRingNumber,
    /// found `[..]` that did not contain an element
    MissingBracketElement,
    /// Missing Element
    MissingElement,
    /// Non Bond in Bracket
    NonBondInBracket,
    /// Ring Number Overflow (greater than 99)
    RingNumberOverflow(u8),
    /// Unexpectedly inside of brackets
    UnexpectedBracketedState,
    /// Unexpected end of string
    UnexpectedEndOfString,
    /// Error indicating that an unexpected character was encountered.
    UnexpectedCharacter(char),
    /// An unexpected `:` has been found
    UnexpectedColon,
    /// An unexpected `-` has been found
    UnexpectedDash,
    /// An unexpected `%` has been found
    UnexpectedPercent,
    /// An unexpected left bracket `[` was found
    UnexpectedLeftBracket,
    /// An unexpected right bracket `]` was found
    UnexpectedRightBracket,
    /// A closing `]` bracket was not found
    UnclosedBracket,
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use SmilesError::{
            BondInBracket, ChargeOverflow, ChargeUnderflow, ElementRequiresBrackets, ElementsRs,
            IntegerOverflow, InvalidAromaticElement, InvalidChirality, InvalidClass,
            InvalidElementName, InvalidIsotope, InvalidNumber, InvalidRingNumber,
            InvalidUnbracketedAtom, MissingBracketElement, MissingElement, NonBondInBracket,
            RingNumberOverflow, UnclosedBracket, UnexpectedBracketedState, UnexpectedCharacter,
            UnexpectedColon, UnexpectedDash, UnexpectedEndOfString, UnexpectedLeftBracket,
            UnexpectedPercent, UnexpectedRightBracket,
        };
        match self {
            MissingElement => write!(f, "Missing element"),
            InvalidIsotope => write!(f, "Invalid isotope"),
            InvalidElementName(c) => write!(f, "Invalid element name: {c}"),
            InvalidNumber => write!(f, "Invalid number"),
            UnexpectedCharacter(c) => write!(f, "Unexpected character: {c}"),
            UnexpectedLeftBracket => write!(f, "Unexpected '['"),
            UnexpectedRightBracket => write!(f, "Unexpected ']'"),
            UnclosedBracket => write!(f, "Unclosed '['"),
            ElementRequiresBrackets => write!(f, "Element requires brackets"),
            MissingBracketElement => write!(f, "Missing element inside brackets"),
            InvalidAromaticElement(e) => write!(f, "Invalid aromatic element: {e}"),
            IntegerOverflow => write!(f, "Integer overflow"),
            RingNumberOverflow(n) => write!(f, "Ring number overflow: {n}"),
            ChargeUnderflow(c) => write!(f, "Charge underflow: {c}"),
            ChargeOverflow(c) => write!(f, "Charge overflow: {c}"),
            BondInBracket(b) => write!(f, "Bond in bracket: {b}"),
            NonBondInBracket => write!(f, "Non-bond '.' in bracket"),
            InvalidChirality => write!(f, "Invalid chirality"),
            UnexpectedEndOfString => write!(f, "Unexpected end of string"),
            InvalidClass => write!(f, "Invalid class"),
            InvalidUnbracketedAtom(a) => write!(f, "Invalid unbracketed atom: {a}"),
            UnexpectedBracketedState => write!(f, "Unexpected bracketed state"),
            UnexpectedDash => write!(f, "Unexpected '-'"),
            UnexpectedColon => write!(f, "Unexpected ':'"),
            UnexpectedPercent => write!(f, "Unexpected '%'"),
            InvalidRingNumber => write!(f, "Invalid ring number"),
            ElementsRs(error) => write!(f, "Error Parsing Element: {error}"),
        }
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
#[derive(Debug)]
pub struct SmilesErrorWithSpan {
    /// The [`SmilesError`]
    smiles_error: SmilesError,
    /// The span as `usize`
    span: Range<usize>,
}

impl SmilesErrorWithSpan {
    /// Creates a new error from the [`SmilesError`] and the `span`
    #[must_use]
    pub fn new(smiles_error: SmilesError, start: usize, end: usize) -> Self {
        Self { smiles_error, span: Range { start, end } }
    }

    /// Returns the [`SmilesError`]
    #[must_use]
    pub fn smiles_error(&self) -> SmilesError {
        self.smiles_error
    }

    /// Returns the start of the span
    #[must_use]
    pub fn start(&self) -> usize {
        self.span.start
    }

    /// Returns the end of the span
    #[must_use]
    pub fn end(&self) -> usize {
        self.span.end
    }

    /// Returns the full span for the error
    #[must_use]
    pub fn span(&self) -> &Range<usize> {
        &self.span
    }

    /// Render the error pointing back to location in the original string
    #[must_use]
    pub fn render(&self, input: &str) -> String {
        let start = self.start().min(input.len());
        let end = self.end().min(input.len()).max(start + 1).min(input.len());

        let mut underline = String::new();
        underline.push_str(&" ".repeat(start));
        underline.push_str(&"^".repeat(end - start));

        format!("{input}\n{underline}\n{}", self.smiles_error)
    }
}

impl fmt::Display for SmilesErrorWithSpan {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} at {}..{}", self.smiles_error, self.start(), self.end())
    }
}
