# SMILES-parser
[![Rust CI](https://github.com/earth-metabolome-initiative/smiles-parser/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/earth-metabolome-initiative/smiles-parser/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/earth-metabolome-initiative/smiles-parser/graph/badge.svg)](https://codecov.io/gh/earth-metabolome-initiative/smiles-parser)

A parser that checks the validity of SMILES strings and converts them into molecular graph representations.

## Parsing Specification
This parser was designed by following the [OpenSMILES specification](http://opensmiles.org/opensmiles.html) and [Wikipedia Article](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System).

## SMILES Parsing Rules:

### Valid SMILES Characters:

| Character | SMILES Purpose |
| --------- | -------------- |
| `A–Z, a–z` | **Atom symbols**. Unbracketed atoms are limited to the organic subset `B C N O P S F Cl Br I` plus aromatic forms `b c n o p s` and `*`. Inside brackets, atom symbols must match a valid element symbol (e.g., `Si`, `Na`, `se`, `as`) or `*`. Aromatic atoms are denoted by lowercase symbols. Multi-character element symbols use initial uppercase + lowercase (e.g., `Cl`, `Br`, `Si`), except bracketed aromatic symbols like `se`, `as`. | 
| `[` `]` | **Bracket atoms**. Enter/exit bracket-atom grammar: `[` *isotope? symbol chiral? hcount? charge? class?* `]`. Brackets are required for non-organic-subset elements and whenever isotope, explicit hydrogen count, charge, chirality, or atom class are specified. |
| `*` | **Wildcard atom**. May appear unbracketed (`*`) or bracketed (`[*]`), and in brackets may carry isotope/chirality/H-count/charge/class. |
| `@` | **Chirality tag introducer** inside bracket atoms. Used as `@` / `@@` and extended forms like `@TH1`, `@AL1`, `@SP1`, `@TB1`, `@OH1` |
| `+` `-` | **Charge signs** inside bracket atoms (e.g., `[O-]`, `[Cu+2]`, `[Ti++++]`). Note: `-` is also a **bond symbol** outside brackets.|
| `:` `-` `=` `#` `$` `/` `\` `.` | **Bond symbols** in the main chain. `.` is the dot/disconnect (“no bond between components”). `-` is an explicit single bond (assumed to be single bond if bond is omitted from notation) and must be distinguished from charge sign by context (inside vs outside brackets). `:` represents an aromatic *one and a half* bond but may also be used for class. |
| `:` `0-9` | The `:` may also be used to represent arbitrary integers that do not have chemical meaning in the SMILES string (inside of brackets only), but may be used by applications working with SMILES strings, for classifying atoms in said applications (`[CH4:2]` marks Methane as being in class `2`)j.|  
| `%` `0–9` | **Digits** occur in multiple sub-grammars. Outside brackets, digits denote **ring closures**: `0–9` for single-digit ring numbers, and `%` followed by **exactly two digits** for ring numbers `00–99` (e.g., `C%12...%12`). Inside brackets, digits may appear as **isotope** (before symbol), **H-count** (after `H`), **charge magnitude** (after `+`/`-`), and **atom class** (after `:`). Note: `%123` is parsed as ring closure `%12` followed by ring closure `3`. |
| `(` `)` | **Branching**. Parentheses introduce a branch off the current atom. |

## Example Usage

```rust
use core::str::FromStr;

use smiles_parser::smiles::Smiles;

let smiles = Smiles::from_str("CCO").expect("valid SMILES should parse");

// The molecular graph contains 3 atoms and 2 bonds for ethanol.
assert_eq!(smiles.nodes().len(), 3);
assert_eq!(smiles.number_of_bonds(), 2);

// You can inspect each node in the graph.
let node_summaries: Vec<(usize, String)> = smiles
    .nodes()
    .iter()
    .enumerate()
    .map(|(index, atom)| (index, atom.to_string()))
    .collect();

assert_eq!(
    node_summaries,
    vec![
        (0, "C".to_string()),
        (1, "C".to_string()),
        (2, "O".to_string()),
    ]
);

// You can also inspect the rendered graph again as a SMILES string.
assert_eq!(smiles.render().unwrap(), "CCO");
```
