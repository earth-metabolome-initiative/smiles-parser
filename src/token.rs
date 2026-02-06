//! Represents tokens used in parsing SMILES strings.

use elements_rs::Element;

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// An unknown/wild card atom
    Asterisk,
    /// An atom inside brackets `[]`s
    Bracketed(BracketedAtom),
    /// A bond between atoms
    Bond(Bond),
    /// A full stop `.`
    Dot,
    /// A left parentheses `(`
    LeftParentheses,
    /// An percent sign `%`
    Percent,
    /// A right parentheses `)`
    RightParentheses,
    /// An atom occurring outside of brackets `[]`
    Unbracketed(UnbracketedAtom)
}

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct UnbracketedAtom {
    /// Unbracketed elements as [`Element`]
    pub element: Element,
    /// Whether the atom is aromatic
    pub aromatic: bool,
}

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct BracketedAtom {
    /// Bracketed elements as [`Element`]
    pub element: Element,
    /// Parsed Isotope Mass Number Value
    pub isotope_mass_number: Option<u16>,
    pub aromatic: bool,
    pub hydrogens: HydrogenCount,
    pub charge: Charge,
}

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub enum HydrogenCount {
    Unspecified,
    Explicit(u8),
}

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Charge(pub Option<i8>);

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub enum Bond {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Up,
    Down,
}