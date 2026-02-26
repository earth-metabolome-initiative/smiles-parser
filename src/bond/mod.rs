//! Module for specifying the bond between two atoms in a `SMILES` string
pub mod bond_edge;
pub mod ring_num;


use core::fmt;

#[derive(Copy, Debug, Default, PartialEq, Clone, Eq, Hash)]
/// Enum used to specify the Bond type, based on SMILES specification
pub enum Bond {
    #[default]
    /// Implicit single bond or explicit with `-`
    Single,
    /// Defined with `=`
    Double,
    /// Defined with `#`
    Triple,
    /// Defined with `$`
    Quadruple,
    /// Aromatic bonds defined with `:`
    Aromatic,
    /// Represents a stereochemical single bond `/` (up)
    Up,
    /// Represents a stereochemical single bond `\` (down)
    Down,
}

impl fmt::Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Bond::Single => "-",
            Bond::Double => "=",
            Bond::Triple => "#",
            Bond::Quadruple => "$",
            Bond::Aromatic => ":",
            Bond::Up => "/",
            Bond::Down => "\\",
        };
        f.write_str(s)
    }
}

#[cfg(test)]
mod tests {
    use crate::bond::Bond;

    #[test]
    fn test_default() {
        assert_eq!(Bond::default(), Bond::Single);
    }
}
