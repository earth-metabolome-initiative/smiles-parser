#![no_main]

use libfuzzer_sys::fuzz_target;
use smiles_parser::smiles::Smiles;

fuzz_target!(|data: &str| {
    if let Ok(smiles) = data.parse::<Smiles>() {
        let rendered = smiles.to_string();

        if let Ok(reparsed) = rendered.parse::<Smiles>() {
            let rerendered = reparsed.to_string();
            assert_eq!(&rendered, &rerendered,);
        }
    }
});

