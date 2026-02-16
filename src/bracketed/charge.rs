//! Module for validating a charge on an atom
use crate::errors::SmilesError;

#[derive(Copy, Debug, PartialEq, Clone, Eq, Hash)]
/// Wrapper struct for possible charge to limit upper and lower bounds
pub struct Charge(i8);

impl Charge {
    /// Attempts to set the `Charge`, if outside of bounds returns
    /// [`SmilesError::ChargeUnderflow`] or [`SmilesError::ChargeOverflow`]
    pub fn try_new(num: i8) -> Result<Self, SmilesError> {
        (-15..=15).contains(&num).then_some(Self(num)).ok_or(match num.is_negative() {
            true => SmilesError::ChargeUnderflow(num),
            false => SmilesError::ChargeOverflow(num),
        })
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
