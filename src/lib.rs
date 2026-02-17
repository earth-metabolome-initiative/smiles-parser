#![doc = include_str!("../README.md")]

pub mod atom_symbol;
pub mod bond;
pub mod bracketed;
pub mod errors;
pub mod parser;
pub mod ring_num;
pub mod smiles;
pub mod token;
pub mod unbracketed;
/// A prelude module to simplify imports.
pub mod prelude {
    pub use crate::smiles::Smiles;
}
