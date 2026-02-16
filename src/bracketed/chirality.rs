//! Module for parsing, validating, and specifying the chirality of an atom
use crate::errors::SmilesError;

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

#[cfg(test)]
mod tests {
    use crate::bracketed::chirality::Chirality;

    #[test]
    fn test_chiral_range_bounds() {
        let th_err = Chirality::try_th(3);
        assert_eq!(th_err, Err(crate::errors::SmilesError::InvalidChirality));
    }
}
