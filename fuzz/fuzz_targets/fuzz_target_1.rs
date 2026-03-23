#![no_main]

use libfuzzer_sys::fuzz_target;
use std::str;
use smiles_parser::smiles::Smiles;

fuzz_target!(|data: &str| {
    if let Ok(smiles) = data.parse::<Smiles>() {
        let rendered = smiles.to_string();
        if let Ok(smiles_reparsed) = rendered.parse::<Smiles>() {
        assert_eq!(smiles, smiles_reparsed);
        }
    }
});
