//! Module for specifying the total number of hydrogens a `SMILES` string
//! specifies

use std::fmt;

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
    #[must_use]
    pub fn new(hydrogens: Option<u8>) -> Self {
        match hydrogens {
            Some(h) => Self::Explicit(h),
            None => Self::Unspecified,
        }
    }
    /// Retrieves the hydrogen count if explicit
    #[must_use]
    pub fn get_count(&self) -> Option<u8> {
        match self {
            Self::Unspecified => None,
            Self::Explicit(n) => Some(*n),
        }
    }
}

impl fmt::Display for HydrogenCount {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unspecified => Ok(()),
            Self::Explicit(n) => write!(f, "H{n}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::atom::bracketed::hydrogen_count::HydrogenCount;

    #[test]
    fn test_hydro_variants() {
        let count = HydrogenCount::new(Some(10));
        assert_eq!(count.get_count(), Some(10));

        let no_count = HydrogenCount::default();
        assert_eq!(no_count.get_count(), None);

        let no_count = HydrogenCount::new(None);
        assert_eq!(no_count.get_count(), None);
        assert_eq!(no_count, HydrogenCount::Unspecified);
        assert_eq!(no_count, HydrogenCount::default());
    }

    #[test]
    fn test_hydrogen_count_fmt_all_arms() {
        let no_count = HydrogenCount::Unspecified;
        let count = HydrogenCount::Explicit(5);
        assert_eq!(no_count.to_string(), String::new());
        assert_eq!(count.to_string(), "H5");
    }
}
