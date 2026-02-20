//! Module for validating a charge on an atom
use crate::errors::SmilesError;

#[derive(Copy, Default, Debug, PartialEq, Clone, Eq, Hash)]
/// Wrapper struct for possible charge to limit upper and lower bounds
pub struct Charge(i8);

impl Charge {
    /// Attempts to set the `Charge`, if outside of bounds returns
    /// [`SmilesError::ChargeUnderflow`] or [`SmilesError::ChargeOverflow`]
    ///
    /// # Errors
    /// - Returns [`SmilesError::ChargeUnderflow`] if `i8` is less than `-15`
    /// - Returns [`SmilesError::ChargeOverflow`] if `i8` is greater than `15`
    pub fn try_new(num: i8) -> Result<Self, SmilesError> {
        (-15..=15).contains(&num).then_some(Self(num)).ok_or(if num.is_negative() {
            SmilesError::ChargeUnderflow(num)
        } else {
            SmilesError::ChargeOverflow(num)
        })
    }

    /// Returns the `Charge` value as `i8`
    #[must_use]
    pub fn get(&self) -> i8 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::Charge;
    use crate::errors::SmilesError;

    #[test]
    fn default_is_zero() {
        assert_eq!(Charge::default().get(), 0);
    }

    #[test]
    fn try_new_accepts_in_range_values() {
        assert_eq!(Charge::try_new(-15).map(|c| c.get()), Ok(-15));
        assert_eq!(Charge::try_new(15).map(|c| c.get()), Ok(15));
        assert_eq!(Charge::try_new(-1).map(|c| c.get()), Ok(-1));
        assert_eq!(Charge::try_new(0).map(|c| c.get()), Ok(0));
        assert_eq!(Charge::try_new(1).map(|c| c.get()), Ok(1));
        assert_eq!(Charge::try_new(7).map(|c| c.get()), Ok(7));
        assert_eq!(Charge::try_new(-7).map(|c| c.get()), Ok(-7));
    }

    #[test]
    fn try_new_rejects_negative_underflow() {
        assert_eq!(Charge::try_new(-16), Err(SmilesError::ChargeUnderflow(-16)));
        assert_eq!(Charge::try_new(i8::MIN), Err(SmilesError::ChargeUnderflow(i8::MIN)));
    }

    #[test]
    fn try_new_rejects_positive_overflow() {
        assert_eq!(Charge::try_new(16), Err(SmilesError::ChargeOverflow(16)));
        assert_eq!(Charge::try_new(i8::MAX), Err(SmilesError::ChargeOverflow(i8::MAX)));
    }

    #[test]
    fn try_new_error_type_depends_on_sign() {
        match Charge::try_new(-16) {
            Err(SmilesError::ChargeUnderflow(-16)) => {}
            other => panic!("expected ChargeUnderflow(-16), got {other:?}"),
        }

        match Charge::try_new(16) {
            Err(SmilesError::ChargeOverflow(16)) => {}
            other => panic!("expected ChargeOverflow(16), got {other:?}"),
        }
    }
}
