#![doc = include_str!("../README.md")]

pub mod atom;
pub mod bond;
pub mod errors;
pub mod parser;
pub mod smiles;
pub mod token;
/// A prelude module to simplify imports.
pub mod prelude {
    pub use crate::smiles::Smiles;
}
