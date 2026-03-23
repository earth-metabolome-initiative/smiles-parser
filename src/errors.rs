//! Defines errors used in the SMILES parser.

use core::fmt;
use std::{num::TryFromIntError, ops::Range};

use elements_rs::Element;

use crate::{atom::atom_symbol::AtomSymbol, bond::Bond};

/// The errors that could occur during SMILES parsing.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SmilesError {
    /// Bond Inside Bracket
    BondInBracket(Bond),
    /// A charge is over the allowed maximum (15)
    ChargeOverflow(i8),
    /// A charge is below allowed minimum (-15)
    ChargeUnderflow(i8),
    /// A Duplicate [`crate::atom::atom_node::AtomNode`] id was found
    DuplicateNodeId(usize),
    /// A non bare element found outside of brackets
    ElementRequiresBrackets,
    /// Wrapper for `element_rs` errors
    ElementsRs(elements_rs::errors::Error),
    /// A bond was not able to bind two atoms
    IncompleteBond(Bond),
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
    /// Invalid `Token::NonBond`
    InvalidNonBondToken,
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
    /// Node id does not point to a valid `AtomNode`
    NodeIdInvalid(usize),
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
    /// A `(` that wasn't expected has been found
    UnexpectedLeftParentheses,
    /// An unexpected right bracket `]` was found
    UnexpectedRightBracket,
    /// An unexpected right parentheses `)` was found
    UnexpectedRightParentheses,
    /// A closing `]` bracket was not found
    UnclosedBracket,
    /// A branch has not been closed with a `)`
    UnclosedBranch,
    /// A ring number has been found that was not completed
    UnclosedRing,
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use SmilesError::{
            BondInBracket, ChargeOverflow, ChargeUnderflow, DuplicateNodeId,
            ElementRequiresBrackets, ElementsRs, IncompleteBond, IntegerOverflow,
            InvalidAromaticElement, InvalidChirality, InvalidClass, InvalidElementName,
            InvalidIsotope, InvalidNonBondToken, InvalidNumber, InvalidRingNumber,
            InvalidUnbracketedAtom, MissingBracketElement, MissingElement, NodeIdInvalid,
            NonBondInBracket, RingNumberOverflow, UnclosedBracket, UnclosedBranch, UnclosedRing,
            UnexpectedBracketedState, UnexpectedCharacter, UnexpectedColon, UnexpectedDash,
            UnexpectedEndOfString, UnexpectedLeftBracket, UnexpectedLeftParentheses,
            UnexpectedPercent, UnexpectedRightBracket, UnexpectedRightParentheses,
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
            InvalidNonBondToken => write!(f, "Invalid Non-bond '.' found"),
            NodeIdInvalid(n) => write!(f, "Invalid Atom Node ID: {n}"),
            IncompleteBond(bond) => write!(f, "Bond: {bond} missing Atom Node(s)"),
            UnexpectedLeftParentheses => write!(f, "Unexpected '('"),
            UnclosedBranch => write!(f, "Branch not closed"),
            UnclosedRing => write!(f, "Ring not closed"),
            UnexpectedRightParentheses => write!(f, "Unexpected `)`"),
            DuplicateNodeId(id) => write!(f, "Node ID: {id} duplicated"),
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
    pub fn span(&self) -> Range<usize> {
        self.span.start..self.span.end
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

#[cfg(test)]
mod tests {
    use std::num::TryFromIntError;

    use elements_rs::Element;

    use crate::{
        atom::atom_symbol::AtomSymbol,
        bond::Bond,
        errors::{SmilesError, SmilesErrorWithSpan},
    };

    #[test]
    fn test_smiles_error_display_monolithic() {
        let elements_rs_error = elements_rs::errors::Error::AtomicNumber(4);

        let cases = [
            (
                SmilesError::BondInBracket(Bond::Aromatic),
                format!("Bond in bracket: {}", Bond::Aromatic),
            ),
            (SmilesError::ChargeOverflow(50), "Charge overflow: 50".to_string()),
            (SmilesError::ChargeUnderflow(-50), "Charge underflow: -50".to_string()),
            (SmilesError::DuplicateNodeId(2), "Node ID: 2 duplicated".to_string()),
            (SmilesError::ElementRequiresBrackets, "Element requires brackets".to_string()),
            (
                SmilesError::ElementsRs(elements_rs_error),
                format!("Error Parsing Element: {elements_rs_error}"),
            ),
            (
                SmilesError::IncompleteBond(Bond::Aromatic),
                format!("Bond: {} missing Atom Node(s)", Bond::Aromatic),
            ),
            (
                SmilesError::InvalidAromaticElement(Element::Ac),
                format!("Invalid aromatic element: {}", Element::Ac),
            ),
            (SmilesError::InvalidChirality, "Invalid chirality".to_string()),
            (SmilesError::InvalidClass, "Invalid class".to_string()),
            (SmilesError::InvalidElementName('w'), "Invalid element name: w".to_string()),
            (SmilesError::InvalidIsotope, "Invalid isotope".to_string()),
            (SmilesError::InvalidNonBondToken, "Invalid Non-bond '.' found".to_string()),
            (SmilesError::InvalidNumber, "Invalid number".to_string()),
            (SmilesError::IntegerOverflow, "Integer overflow".to_string()),
            (
                SmilesError::InvalidUnbracketedAtom(AtomSymbol::WildCard),
                format!("Invalid unbracketed atom: {}", AtomSymbol::WildCard),
            ),
            (SmilesError::InvalidRingNumber, "Invalid ring number".to_string()),
            (SmilesError::MissingBracketElement, "Missing element inside brackets".to_string()),
            (SmilesError::MissingElement, "Missing element".to_string()),
            (SmilesError::NodeIdInvalid(2), "Invalid Atom Node ID: 2".to_string()),
            (SmilesError::NonBondInBracket, "Non-bond '.' in bracket".to_string()),
            (SmilesError::RingNumberOverflow(100), "Ring number overflow: 100".to_string()),
            (SmilesError::UnexpectedBracketedState, "Unexpected bracketed state".to_string()),
            (SmilesError::UnexpectedEndOfString, "Unexpected end of string".to_string()),
            (SmilesError::UnexpectedCharacter('$'), "Unexpected character: $".to_string()),
            (SmilesError::UnexpectedColon, "Unexpected ':'".to_string()),
            (SmilesError::UnexpectedDash, "Unexpected '-'".to_string()),
            (SmilesError::UnexpectedPercent, "Unexpected '%'".to_string()),
            (SmilesError::UnexpectedLeftBracket, "Unexpected '['".to_string()),
            (SmilesError::UnexpectedLeftParentheses, "Unexpected '('".to_string()),
            (SmilesError::UnexpectedRightBracket, "Unexpected ']'".to_string()),
            (SmilesError::UnexpectedRightParentheses, "Unexpected `)`".to_string()),
            (SmilesError::UnclosedBracket, "Unclosed '['".to_string()),
            (SmilesError::UnclosedBranch, "Branch not closed".to_string()),
            (SmilesError::UnclosedRing, "Ring not closed".to_string()),
        ];

        for (error, expected) in cases {
            assert_eq!(error.to_string(), expected, "wrong Display for {error:?}");
        }
    }

    #[test]
    fn test_smiles_error_from_conversions() {
        let int_error: Result<u8, TryFromIntError> = u8::try_from(300_i16);
        let int_error = int_error.unwrap_err();

        let smiles_error: SmilesError = int_error.into();
        assert_eq!(smiles_error, SmilesError::InvalidNumber);

        let elements_error = elements_rs::errors::Error::AtomicNumber(4);
        let smiles_error: SmilesError = elements_error.into();
        assert_eq!(
            smiles_error,
            SmilesError::ElementsRs(elements_rs::errors::Error::AtomicNumber(4))
        );
    }

    #[test]
    fn test_smiles_error_with_span() {
        let error = SmilesErrorWithSpan::new(SmilesError::UnexpectedCharacter('$'), 2, 3);

        assert_eq!(error.smiles_error(), SmilesError::UnexpectedCharacter('$'));
        assert_eq!(error.start(), 2);
        assert_eq!(error.end(), 3);
        assert_eq!(error.span(), (2..3));

        assert_eq!(error.to_string(), "Unexpected character: $ at 2..3");

        assert_eq!(error.render("CC$O"), "CC$O\n  ^\nUnexpected character: $");
    }
}
