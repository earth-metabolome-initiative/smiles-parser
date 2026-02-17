//! Module for mapping and validating a ring marker
use crate::errors::SmilesError;


#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Represents a ring marker and implements tighter bounds for the minimal and
/// maximal value a ring marker can be
pub struct RingNum(u8);
impl RingNum {
    /// Attempts to generate a [`RingNum`] form a [`u8`],
    ///
    /// # Errors
    /// - Returns a [`SmilesError::RingNumberOverflow`] if the value is above
    ///   `99`
    pub fn try_new(num: u8) -> Result<Self, SmilesError> {
        (0..=99).contains(&num).then_some(Self(num)).ok_or(SmilesError::RingNumberOverflow(num))
    }

    /// Returns the value for the [`RingNum`]
    pub fn get(&self) -> u8 {
        self.0
    }
}
