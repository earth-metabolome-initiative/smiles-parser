use std::str::FromStr;

use super::Smiles;
use crate::{errors::SmilesErrorWithSpan, parser::token_iter::TokenIter};

impl FromStr for Smiles {
    type Err = SmilesErrorWithSpan;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let token_iter = TokenIter::from(s);
        let _tokens = token_iter.collect::<Result<Vec<_>, _>>()?;
        todo!()
    }
}
