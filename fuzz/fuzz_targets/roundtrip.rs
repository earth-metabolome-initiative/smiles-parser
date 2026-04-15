#![no_main]

use std::string::ToString;

use libfuzzer_sys::fuzz_target;
use smiles_parser::smiles::Smiles;


fuzz_target!(|data: &str| {
    match data.parse::<Smiles>() {
        Ok(smiles) => {
            let rendered = smiles.to_string();
            let reparsed = rendered.parse::<Smiles>().unwrap();
            let rerendered = reparsed.to_string();
            assert_eq!(
                rendered, rerendered,
                "rendered forms differ after round-trip: rendered={rendered} rerendered={rerendered}"
            );
        }
        Err(err) => {
            let _ = err.render(data);
        }
    }
});
