//! Represents tokens used in parsing SMILES strings.

use elements_rs::{Element, Isotope};

#[derive(Debug, PartialEq, Clone, Eq, PartialOrd, Ord, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// the ampersand `&`
    Ampersand,
    /// An Atom with associated properties
    Atom(AtomToken),
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
    /// A left parentheses `(`
    LeftParentheses,
    /// The minus sign `-`
    Minus, 
    /// The percent sign `%`
    Percent, 
    /// The plus sign `+`
    Plus,
    /// A right parentheses `)`
    RightParentheses,
}

#[derive(Debug, PartialEq, Clone, Eq, PartialOrd, Ord, Hash)]
/// Represents the Atomic Token variants in a SMILES string
pub enum AtomToken {
    /// The atoms in the [organic subset](http://opensmiles.org/opensmiles.html#orgsbst)
    Bare {
        /// The element as an [`Element`] 
        element: Element,
        /// Whether the bare atom is aromatic   
        aromatic: bool,     
    },
    /// Atoms parsed from within square brackets `[..]` 
    Bracketed {
        /// Optional isotope specification as [`Isotope`]
        isotope: Option<Isotope>,
        /// The element as an [`Element`]
        element: Element,
        /// Whether the atom is aromatic
        aromatic: bool,
        /// The charge of the atom as an [`i8`]
        charge: i8,
        /// The number of bonded hydrogens as [`u8`]
        hydrogens: u8,
    },
}
