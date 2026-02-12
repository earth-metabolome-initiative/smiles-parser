//! Represents tokens used in parsing SMILES strings.

use std::{default, ops::Range};

use elements_rs::{Element, Isotope};

use crate::errors::SmilesError;

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a token in a molecular formula.
pub enum Token {
    /// Represented with a `.`
    NonBond,
    /// Elements that occur inside of `[]`, structured as [`BracketedAtom`]
    BracketedAtom(BracketedAtom),
    /// Aliphatic organic elements only as [`UnbracketedAtom`]
    UnbracketedAtom(UnbracketedAtom),
    /// The parsed bond. Single bonds are implicit between atoms if not
    /// explicated as a bond, e.g. `CC` would be `C-C`.
    Bond(Bond),
    /// Token for left parentheses, used for branching `(`
    LeftParentheses,
    /// Token for right parentheses, used for ending a branch `)`
    RightParentheses,
    /// Ring number markers occur outside of `[]` and may be of type `%` and
    /// `0-99`, and `%` may be omitted.
    RingClosure(RingNum),
}

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
/// A parsed token and its relative location in the string
pub struct TokenWithSpan {
    /// The parsed token
    token: Token,
    /// The startin to ending points of the token in the string
    span: Range<usize>,
}

impl TokenWithSpan {
    /// Generates a new token with a specified span
    pub fn new(token: Token, start: usize, end: usize) -> Self {
        Self { token, span: Range { start, end } }
    }
    /// Returns the token
    pub fn token(&self) -> Token {
        self.token
    }
    /// Returns the span as [`Range`]
    pub fn span(&self) -> &Range<usize> {
        &self.span
    }
    /// Returns the start of the span as [`usize`]
    pub fn start(&self) -> usize {
        self.span.start
    }
    /// Returns the end of the span as [`usize`]
    pub fn end(&self) -> usize {
        self.span.end
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a ring marker and implements tighter bounds for the minimal and
/// maximal value a ring marker can be
pub struct RingNum(u8);
impl RingNum {
    /// Maximum value for a ring marker
    pub const MAX: u8 = 99;

    /// Attempts to generate a [`RingNum`] form a [`u8`],
    ///
    /// # Errors
    /// - Returns a [`SmilesError::RingNumberOverflow`] if the value is above
    ///   `99`
    pub fn try_new(n: u8) -> Result<Self, SmilesError> {
        if n <= Self::MAX { Ok(Self(n)) } else { Err(SmilesError::RingNumberOverflow(n)) }
    }

    /// Returns the value for the [`RingNum`]
    pub fn get(&self) -> u8 {
        self.0
    }
}

#[derive(Copy, Default, Debug, PartialEq, Clone, Eq, Hash)]
/// Enum to allow for standard elements or the `WildCard` variant, represented
/// as `*`
pub enum AtomSymbol {
    /// The explicitly named [`Element`]
    Element(Element),
    /// WildCard variant, described [here](http://opensmiles.org/opensmiles.html#inatoms)
    #[default]
    WildCard,
}

impl AtomSymbol {
    /// Creates an atom symbol
    pub fn new(element_type: Option<Element>) -> Self {
        match element_type {
            Some(element) => AtomSymbol::Element(element),
            None => AtomSymbol::default(),
        }
    }
    /// Verifies whether the symbol present is a wildcard
    pub fn is_wildcard(&self) -> bool {
        matches!(self, AtomSymbol::WildCard)
    }
    /// Returns either the [`Element`] or `None` if wildcard
    pub fn element(&self) -> Option<Element> {
        match self {
            AtomSymbol::Element(e) => Some(*e),
            AtomSymbol::WildCard => None,
        }
    }
    /// Consumes the `AtomSymbol` and returns the [`Element`] or `None` if
    /// `WildCard`
    pub fn into_element(self) -> Option<Element> {
        match self {
            AtomSymbol::Element(e) => Some(e),
            AtomSymbol::WildCard => None,
        }
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Structure for aliphatic atoms, aromatic or non aromatic
pub struct UnbracketedAtom {
    /// Unbracketed elements as [`Element`]
    symbol: AtomSymbol,
    /// Whether the atom is aromatic
    aromatic: bool,
}

impl UnbracketedAtom {
    /// Creates a new `UnbracketedAtom`
    pub const fn new(symbol: AtomSymbol, aromatic: bool) -> Self {
        Self { symbol, aromatic }
    }
    /// Returns the [`AtomSymbol`] of the atom
    pub fn symbol(&self) -> AtomSymbol {
        self.symbol
    }
    /// Returns the [`Element`] or `None` if `WildCard`
    pub fn element(&self) -> Option<Element> {
        self.symbol.element()
    }
    /// Returns true of aromatic
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }
    /// Returns true if `AtomSymbol` is [`AtomSymbol::WildCard`]
    pub fn is_wildcard(&self) -> bool {
        self.symbol.is_wildcard()
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Contains [`Element`] and specified meta data about an element in `[]`
pub struct BracketedAtom {
    /// Bracketed elements as [`Element`]
    symbol: AtomSymbol,
    /// Parsed Isotope Mass Number Value
    isotope_mass_number: Option<u16>,
    /// If bracketed element is aromatic
    aromatic: bool,
    /// The number of Hydrogens explicitly listed or `Unspecified`
    hydrogens: HydrogenCount,
    /// The charge of the Atom, default is `0`
    charge: Charge,
    /// Esoteric potential integers presented after Element name: `[CH4:2]`.
    /// Unspecified default is 0
    class: u16,
    /// Denotes Chirality if present
    chiral: Option<Chirality>,
}

impl BracketedAtom {
    /// Returns a builder for bracketed atom
    pub fn builder() -> BracketedAtomBuilder {
        BracketedAtomBuilder {
            bracket_atom: Self {
                symbol: AtomSymbol::default(),
                aromatic: false,
                isotope_mass_number: None,
                hydrogens: HydrogenCount::Unspecified,
                charge: Charge::default(),
                class: 0,
                chiral: None,
            },
        }
    }
    /// Returns the the [`Element`] of the bracket atom or `None` for `WildCard`
    pub fn element(&self) -> Option<Element> {
        match self.symbol {
            AtomSymbol::WildCard => None,
            AtomSymbol::Element(element) => Some(element),
        }
    }
    /// Returns the [`AtomSymbol`]
    pub fn symbol(&self) -> AtomSymbol {
        self.symbol
    }
    /// Returns the isotope mass number
    pub fn isotope_mass_number(&self) -> Option<u16> {
        self.isotope_mass_number
    }
    /// Returns the [`Isotope`] for the [`Element`] for the atom
    pub fn isotope(&self) -> Result<Isotope, SmilesError> {
        let element = self.element().ok_or(SmilesError::InvalidIsotope)?;
        let isotope = match self.isotope_mass_number() {
            None => element.most_abundant_isotope(),
            Some(mass) => Isotope::try_from((element, mass))?,
        };
        Ok(isotope)
    }
    /// Returns aromatic status
    pub fn aromatic(&self) -> bool {
        self.aromatic
    }
    /// Returns the [`HydrogenCount`]
    pub fn hydrogens(&self) -> HydrogenCount {
        self.hydrogens
    }
    /// Returns the hydrogens attached or `None`
    pub fn hydrogen_count(&self) -> Option<u8> {
        match self.hydrogens {
            HydrogenCount::Unspecified => None,
            HydrogenCount::Explicit(i) => Some(i),
        }
    }
    /// Returns the [`Charge`] of the atom
    pub fn charge(&self) -> Charge {
        self.charge
    }
    /// Returns the charge as `i8`
    pub fn charge_value(&self) -> i8 {
        self.charge.get()
    }
    /// Returns the class (default is 0)
    pub fn class(&self) -> u16 {
        self.class
    }
    /// Returns the [`Chirality`] of the atom
    pub fn chiral(&self) -> Option<Chirality> {
        self.chiral
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
/// Builder structure, contains a [`BracketedAtom`] that is mutable until
/// calling `build()`.
pub struct BracketedAtomBuilder {
    bracket_atom: BracketedAtom,
}

impl BracketedAtomBuilder {
    /// Adds an isotope value
    pub fn with_isotope(mut self, iso: u16) -> Self {
        self.bracket_atom.isotope_mass_number = Some(iso);
        self
    }
    /// Adds the element
    pub fn with_element(mut self, possible_element: Option<Element>) -> Self {
        self.bracket_atom.symbol = AtomSymbol::new(possible_element);
        self
    }
    /// Adds the aromatic value
    pub fn with_aromatic(mut self, aromatic: bool) -> Self {
        self.bracket_atom.aromatic = aromatic;
        self
    }
    /// Adds a specified [`HydrogenCount`]
    pub fn with_hydrogens(mut self, h_count: HydrogenCount) -> Self {
        self.bracket_atom.hydrogens = h_count;
        self
    }
    /// Adds a specified [`Charge`]
    pub fn with_charge(mut self, charge: Charge) -> Self {
        self.bracket_atom.charge = charge;
        self
    }
    /// Adds a designated esoteric class
    pub fn with_class(mut self, class: u16) -> Self {
        self.bracket_atom.class = class;
        self
    }
    /// Adds a specified [`Chirality`]
    pub fn with_chiral(mut self, chiral: Chirality) -> Self {
        self.bracket_atom.chiral = Some(chiral);
        self
    }
    /// Consumes the builder and returns the completed [`BracketedAtom`]
    pub fn build(self) -> BracketedAtom {
        BracketedAtom {
            symbol: self.bracket_atom.symbol,
            aromatic: self.bracket_atom.aromatic,
            isotope_mass_number: self.bracket_atom.isotope_mass_number,
            hydrogens: self.bracket_atom.hydrogens,
            charge: self.bracket_atom.charge,
            class: self.bracket_atom.class,
            chiral: self.bracket_atom.chiral,
        }
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Specifies the chirality if present
pub enum Chirality {
    /// `@`
    At,
    /// `@@`
    AtAt,
    /// `@TH` variants (1-2)
    TH(u8),
    /// `@AL` variants (1-2)
    AL(u8),
    /// `@SP` variants (1-3)
    SP(u8),
    /// `@TB` variants (1-20)
    TB(u8),
    /// `@OH` variants (1-30)
    OH(u8),
}

impl Chirality {
    /// Convert `u8` to `TH`+`U8`
    pub fn try_th(num: u8) -> Result<Self, SmilesError> {
        (1..=2).contains(&num).then_some(Self::TH(num)).ok_or(SmilesError::InvalidChirality)
    }
    /// Convert `u8` to `AL`+`U8
    pub fn try_al(num: u8) -> Result<Self, SmilesError> {
        (1..=2).contains(&num).then_some(Self::AL(num)).ok_or(SmilesError::InvalidChirality)
    }
    /// Convert `u8` to `SP`+`U8
    pub fn try_sp(num: u8) -> Result<Self, SmilesError> {
        (1..=3).contains(&num).then_some(Self::SP(num)).ok_or(SmilesError::InvalidChirality)
    }
    /// Convert `u8` to `TB`+`U8
    pub fn try_tb(num: u8) -> Result<Self, SmilesError> {
        (1..=20).contains(&num).then_some(Self::TB(num)).ok_or(SmilesError::InvalidChirality)
    }
    /// Convert `u8` to `OH`+`U8
    pub fn try_oh(num: u8) -> Result<Self, SmilesError> {
        (1..=30).contains(&num).then_some(Self::OH(num)).ok_or(SmilesError::InvalidChirality)
    }
}

#[derive(Copy, Default, Debug, PartialEq, Clone, Eq, Hash)]
/// Designates the hydrogen count (explicit only). Currently Hydrogen count has
/// no upper bound, and may go to [`u8::MAX`]
pub enum HydrogenCount {
    /// Defaults to unspecified
    #[default]
    Unspecified,
    /// The explicit number of hydrogens
    Explicit(u8),
}

impl HydrogenCount {
    /// Returns the hydrogen count based on `Option` input of `u8`
    pub fn new(hydrogens: Option<u8>) -> Self {
        match hydrogens {
            Some(h) => Self::Explicit(h),
            None => Self::Unspecified,
        }
    }
}

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Wrapper struct for possible charge to limit upper and lower bounds
pub struct Charge(i8);

impl Charge {
    /// SMILES specification has lowest charge of -15
    pub const MIN: i8 = -15;
    /// SMILES specification has highest charge of 15
    pub const MAX: i8 = 15;
    /// Attempts to set the `Charge`, if outside of bounds returns
    /// [`SmilesError::ChargeUnderflow`] or [`SmilesError::ChargeOverflow`]
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

    /// Returns the `Charge` value as `i8`
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
/// Enum used to specify the Bond type, based on SMILES specification
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
}

#[cfg(test)]
mod tests {
    use crate::token::Chirality;

    #[test]
    fn test_chiral_range_bounds() {
        let th_err = Chirality::try_th(3);
        assert_eq!(th_err, Err(crate::errors::SmilesError::InvalidChirality));
    }
}
