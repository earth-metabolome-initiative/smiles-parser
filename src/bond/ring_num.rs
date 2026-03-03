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
    #[must_use]
    pub fn get(&self) -> u8 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use crate::{bond::ring_num::RingNum, errors::SmilesError};

    #[test]
    fn test_try_and_get() -> Result<(), SmilesError> {
        let num: u8 = 3;
        let invalid_num: u8 = 200;
        let valid_ring = RingNum::try_new(num)?;
        let invalid_ring = RingNum::try_new(invalid_num);

        assert_eq!(valid_ring.get(), 3);
        assert!(invalid_ring.is_err());

        Ok(())
    }
}
