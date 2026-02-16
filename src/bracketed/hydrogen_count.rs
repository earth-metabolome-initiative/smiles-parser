//! Module for specifying the total number of hydrogens a `SMILES` string
//! specifies
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
