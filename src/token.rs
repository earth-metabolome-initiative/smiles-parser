//! Represents tokens used in parsing SMILES strings.

use std::ops::Range;

use elements_rs::Element;

use crate::errors::SmilesError;

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    BracketedAtom(BracketedAtom),
    UnbracketedAtom(UnbracketedAtom),
    Bond(Bond),
    LeftParentheses,
    RightParentheses,
    RingClosure(RingNum),
}

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct TokenWithSpan {
    token: Token,
    span: Range<usize>,
}

impl TokenWithSpan {
    pub fn new(token: Token, start: usize, end: usize) -> Self {
        Self { token, span: Range { start, end } }
    }
    pub fn token(&self) -> Token {
        self.token
    }
    pub fn span(&self) -> &Range<usize> {
        &self.span
    }
    pub fn start(&self) -> usize {
        self.span.start
    }
    pub fn end(&self) -> usize {
        self.span.end
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
pub struct RingNum(u8);
impl RingNum {
    pub const MIN: u8 = 0;
    pub const MAX: u8 = 99;

    pub fn try_new(n: u8) -> Result<Self, SmilesError> {
        if n <= Self::MAX { Ok(Self(n)) } else { Err(SmilesError::RingNumberOverflow(n)) }
    }

    pub fn get(&self) -> u8 {
        self.0
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
pub enum AtomSymbol {
    Element(Element),
    WildCard,
}

impl AtomSymbol {
    pub fn is_wildcard(&self) -> bool {
        matches!(self, AtomSymbol::WildCard)
    }
    pub fn element(&self) -> Option<&Element> {
        match self {
            AtomSymbol::Element(e) => Some(e),
            AtomSymbol::WildCard => None,
        }
    }
    pub fn into_element(self) -> Option<Element> {
        match self {
            AtomSymbol::Element(e) => Some(e),
            AtomSymbol::WildCard => None,
        }
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
pub struct UnbracketedAtom {
    /// Unbracketed elements as [`Element`]
    symbol: AtomSymbol,
    /// Whether the atom is aromatic
    aromatic: bool,
}

impl UnbracketedAtom {
    pub const fn new(symbol: AtomSymbol, aromatic: bool) -> Self {
        Self { symbol, aromatic }
    }
    pub fn symbol(&self) -> AtomSymbol {
        self.symbol
    }
    pub fn element(&self) -> Option<&Element> {
        self.symbol.element()
    }
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }
    pub fn is_wildcard(&self) -> bool {
        self.symbol.is_wildcard()
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
pub struct BracketedAtom {
    /// Bracketed elements as [`Element`]
    symbol: AtomSymbol,
    /// Parsed Isotope Mass Number Value
    isotope_mass_number: Option<u16>,
    aromatic: bool,
    hydrogens: HydrogenCount,
    charge: Charge,
    /// Esoteric potential integers presented after Element name: `[CH4:2]`.
    /// Unspecified default is 0
    class: u16,
    /// Denotes Chirality if present
    chiral: Option<Chirality>,
}

impl BracketedAtom {
    pub fn builder(symbol: AtomSymbol, aromatic: bool) -> BracketedAtomBuilder {
        BracketedAtomBuilder {
            symbol,
            aromatic,
            isotope_mass_number: None,
            hydrogens: HydrogenCount::Unspecified,
            charge: Charge::default(), 
            class: 0,
            chiral: None,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct BracketedAtomBuilder {
    symbol: AtomSymbol,
    aromatic: bool,
    isotope_mass_number: Option<u16>,
    hydrogens: HydrogenCount,
    charge: Charge,
    class: u16,
    chiral: Option<Chirality>,
}

impl BracketedAtomBuilder {
    pub fn isotope(mut self, iso: u16) -> Self {
        self.isotope_mass_number = Some(iso);
        self
    }

    pub fn hydrogens(mut self, h: HydrogenCount) -> Self {
        self.hydrogens = h;
        self
    }

    pub fn charge(mut self, charge: Charge) -> Self {
        self.charge = charge;
        self
    }

    pub fn class(mut self, class: u16) -> Self {
        self.class = class;
        self
    }

    pub fn chiral(mut self, chiral: Chirality) -> Self {
        self.chiral = Some(chiral);
        self
    }

    pub fn build(self) -> BracketedAtom {
        BracketedAtom {
            symbol: self.symbol,
            aromatic: self.aromatic,
            isotope_mass_number: self.isotope_mass_number,
            hydrogens: self.hydrogens,
            charge: self.charge,
            class: self.class,
            chiral: self.chiral,
        }
    }
}


#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
pub enum Chirality {
    Clockwise,
    AntiClockwise,
}

#[derive(Copy, Default, Debug, PartialEq, Clone, Eq, Hash)]
pub enum HydrogenCount {
    #[default]
    Unspecified,
    Explicit(u8),
}

impl HydrogenCount {
    pub fn new(hydrogens: Option<u8>) -> Self {
        match hydrogens {
            Some(h) => Self::Explicit(h),
            None => Self::Unspecified,
        }
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
pub struct Charge(i8);

impl Charge {
    pub const MIN: i8 = -15;
    pub const MAX: i8 = 15;

    pub fn try_new(n: i8) -> Result<Self, SmilesError> {
        if n >= Self::MIN && n <= Self::MAX {
            Ok(Self(n))
        } else {
            if n.is_negative() {
                Err(SmilesError::ChargeUnderflow(n))
            } else {
                Err(SmilesError::ChargeOverflow(n))
            }
        }
    }
    
    pub fn get(&self) -> i8 {
        self.0
    }
}

impl Default for Charge {
    fn default() -> Self {
        Self(0)
    }
}

#[derive(Copy, Debug, Default, PartialEq, Clone, Eq, Hash)]
pub enum Bond {
    #[default]
    /// Implicit single bond or explicit with `-`
    Single,
    /// Defined with `=`
    Double,
    /// Defined with `#`
    Triple,
    /// Defined with `$`
    Quadruple,
    /// Aromatic bonds defined with `:`
    Aromatic,
    /// Represents a stereochemical single bond `/` (up)
    Up,
    /// Represents a stereochemical single bond `\` (down)
    Down,
    /// Represented with a `.`
    NonBond,
}
