//! Second pass that parses the [`TokenWithSpan`]

use crate::token::TokenWithSpan;

/// Contains the vec of tokens being iterated on and tracks the current position in that vec
pub struct SmilesParser {
    tokens: Vec<TokenWithSpan>,
    position: usize,
}

impl SmilesParser {
    /// Creates a new `SmilesParser` structure
    #[must_use]
    pub fn new(tokens: Vec<TokenWithSpan>) -> Self {
        SmilesParser { tokens, position: 0 }
    }
    /// Retrieves the `tokens` field of [`Vec<TokenWithSpan>`]
    #[must_use]
    pub fn tokens(&self) -> &[TokenWithSpan] {
        &self.tokens
    }
    /// Retrieves the current position
    #[must_use]
    pub fn position(&self) -> usize {
        self.position
    }   
    
    
}
