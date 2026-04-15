#![no_main]

use libfuzzer_sys::fuzz_target;
use smiles_parser::smiles::Smiles;

fuzz_target!(|data: &str| {
    match data.parse::<Smiles>() {
        Ok(smiles) => {
            smiles.debug_assert_canonicalization_invariants();
        }
        Err(err) => {
            let _ = err.render(data);
        }
    }
});
