#![no_main]

use libfuzzer_sys::fuzz_target;
use smiles_parser::smiles::Smiles;

fuzz_target!(|data: &str| {
    if let Ok(smiles) = data.parse::<Smiles>() {
        let rendered = smiles.to_string();

        let reparsed = rendered.parse::<Smiles>().unwrap();
        assert_eq!(smiles, reparsed, "SMILES Graphs are not equivalent smiles 1: {smiles} and smiles 2: {reparsed}, from input: \"{data}\"");
        let rerendered = reparsed.to_string();
        assert_eq!(&rendered, &rerendered, "rendered: {rendered} and rerendered: {rerendered}");
        
    }
});

