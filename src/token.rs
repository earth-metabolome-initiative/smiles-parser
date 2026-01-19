//! Represents tokens used in parsing SMILES strings.

use elements_rs::{Element, Isotope};

#[derive(Debug, PartialEq, Clone, Eq, PartialOrd, Ord, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// An open round bracket
    OpenRoundBracket,
    /// A close round bracket
    CloseRoundBracket,
    /// A dot
    Dot,
    /// A dash '-' character
    Dash,
    /// An equal '=' character i.e. a double bond
    Equal,
    /// A hash '#' character i.e. a triple bond
    Hashtag,
    /// A dollar '$' character i.e. a quadruple bond
    Dollar,
    /// A colon ':' character
    Colon,
    /// A foward slash '/' character
    ForwardSlash,
    /// A back slash '\' character
    BackSlash,
    /// A label that can only go from 0 to 9
    Label(u8),
    ///
    SquareBracketMolecule {
        isotope: Option<Isotope>,
        element: Element,
        aromatic: bool,
        charge: i8,
        hydrogens: u8,
    },
}
