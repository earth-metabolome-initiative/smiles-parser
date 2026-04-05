#![no_std]
#![doc = include_str!("../README.md")]

#[allow(unused_imports)]
#[macro_use]
extern crate alloc;
#[cfg(test)]
#[macro_use]
extern crate std;

pub mod atom;
pub mod bond;
pub mod errors;
pub(crate) mod parser;
pub mod smiles;
pub mod token;
pub(crate) mod traversal;

pub use crate::{
    errors::{SmilesError, SmilesErrorWithSpan},
    smiles::Smiles,
};

/// Common imports for working with this crate.
pub mod prelude {
    pub use crate::{Smiles, SmilesError, SmilesErrorWithSpan};
}
