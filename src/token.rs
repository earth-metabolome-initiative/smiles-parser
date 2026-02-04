//! Represents tokens used in parsing SMILES strings.

use elements_rs::Isotope;

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// the ampersand `&`
    Ampersand,
    /// An unknown/wild card atom
    Asterisk,
    /// A back slash '\' character
    BackSlash,
    /// An atom inside brackets
    BracketAtom(BracketAtom),
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
    /// A left parentheses `(`
    LeftParentheses,
    /// An Organic Atom outside of brackets
    OrganicAtom(Aromatic<OrganicAtom>),
    /// An percent sign `%`
    Percent,
    /// A right parentheses `)`
    RightParentheses,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BracketAtom {
    atoms: Aromatic<Isotope>,
    chirality: Option<Chirality>, 
    hydrogen_count: u8,
    charge: i8,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum OrganicAtom {
    B,
    C,
    N,
    O,
    S,
    P,
    F,
    Cl,
    Br,
    I,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Chirality {
    Clockwise,
    CounterClockwise,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Aromatic<A> {
    atom: A, 
    aromatic: bool, 
}