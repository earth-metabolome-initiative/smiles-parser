//! Represents tokens used in parsing SMILES strings.

use elements_rs::Element;

#[derive(Debug, PartialEq, Clone, Eq, PartialOrd, Ord, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// the ampersand `&`
    Ampersand,
    /// An Atom with associated properties
    Atom {
        /// The Element found
        element: Element,
        /// Whether the element is aromatic
        aromatic: bool,
    },
    /// The at sign `@`
    AtSign,
    /// A back slash '\' character
    BackSlash,
    /// A colon ':' character
    Colon,
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
    /// A left parentheses `(`
    LeftParentheses,
    /// A left square bracket `[`
    LeftSquareBracket,
    /// The minus sign `-`
    Minus,
    /// The percent sign `%`
    Percent,
    /// The plus sign `+`
    Plus,
    /// A right parentheses `)`
    RightParentheses,
    /// A right square bracket `]`
    RightSquareBracket,
}
