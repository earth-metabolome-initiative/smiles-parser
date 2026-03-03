use std::str::FromStr;

use super::Smiles;
use crate::{
    errors::SmilesErrorWithSpan,
    parser::{smiles_parser::SmilesParser, token_iter::TokenIter},
};

impl FromStr for Smiles {
    type Err = SmilesErrorWithSpan;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let token_iter = TokenIter::from(s);
        let tokens = token_iter.collect::<Result<Vec<_>, _>>()?;
        SmilesParser::new(&tokens).parse()
    }
}
