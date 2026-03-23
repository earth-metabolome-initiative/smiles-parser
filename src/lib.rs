#![doc = include_str!("../README.md")]

pub mod atom;
pub mod bond;
pub mod errors;
pub mod parser;
pub mod smiles;
pub mod token;
pub mod traversal;

pub use crate::smiles::Smiles;
pub use crate::errors::{SmilesError, SmilesErrorWithSpan};

/// Common imports for working with this crate.
pub mod prelude {
    pub use crate::{Smiles, SmilesError, SmilesErrorWithSpan};
}