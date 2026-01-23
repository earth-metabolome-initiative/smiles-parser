//! Represents tokens used in parsing SMILES strings.

use elements_rs::{Element, Isotope};

#[derive(Debug, PartialEq, Clone, Eq, PartialOrd, Ord, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// the ampersand `&`
    Ampersand,
    /// The at sign `@`
    AtSign,
    /// A back slash '\' character
    BackSlash,
    /// A colon ':' character
    Colon,
    /// A dash '-' character
    Dash,
    /// A dollar '$' character i.e. a quadruple bond
    Dollar,
    /// A dot
    Dot,
    /// An equal '=' character i.e. a double bond
    Equal,
    /// A forward slash '/' character
    ForwardSlash,
    /// A hash '#' character i.e. a triple bond
    Hashtag,
    /// A label that can only go from 0 to 9
    Label(u8),
    /// A left round bracket `(`
    LeftRoundBracket,
    /// A left square bracket `[`
    LeftSquareBracket,
    /// The minus sign `-`
    Minus, 
    /// The percent sign `%`
    Percent, 
    /// The plus sign `+`
    Plus,
    /// A right round bracket `)`
    RightRoundBracket,
    /// A right square bracket `]`
    RightSquareBracket,
    /// TODO: Figure out if this how we want to shape this enum variant.
    SquareBracketMolecule {
        /// Optional isotope specification
        isotope: Option<Isotope>,
        /// Chemical element symbol inside the brackets.
        element: Element,
        /// Whether the atom is
        aromatic: bool,
        /// Formal charge on the atom
        charge: i8,
        /// Explicit number of hydrogens attached to the atom
        hydrogens: u8,
    },
}
