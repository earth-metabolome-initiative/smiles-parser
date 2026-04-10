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
        f.write_str(self.smiles_symbol())
    }
}

impl Bond {
    #[inline]
    #[must_use]
    pub(crate) const fn smiles_symbol(self) -> &'static str {
        match self {
            Self::Single => "-",
            Self::Double => "=",
            Self::Triple => "#",
            Self::Quadruple => "$",
            Self::Aromatic => ":",
            Self::Up => "/",
            Self::Down => "\\",
        }
    }

    #[inline]
    #[must_use]
    pub(crate) const fn edge_symbol(self) -> &'static str {
        match self {
            Self::Single => "",
            _ => self.smiles_symbol(),
        }
    }

    #[inline]
    #[must_use]
    pub(crate) const fn ring_closure_symbol(self) -> &'static str {
        match self {
            Self::Single | Self::Aromatic => "",
            _ => self.smiles_symbol(),
        }
    }

    #[inline]
    #[must_use]
    pub(crate) const fn without_direction(self) -> Self {
        match self {
            Self::Up | Self::Down => Self::Single,
            _ => self,
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;

    use crate::bond::Bond;

    #[test]
    fn test_default() {
        assert_eq!(Bond::default(), Bond::Single);
    }
    #[test]
    fn test_bond_fmt_all_arms() {
        let cases = [
            (Bond::Single, "-"),
            (Bond::Double, "="),
            (Bond::Triple, "#"),
            (Bond::Quadruple, "$"),
            (Bond::Aromatic, ":"),
            (Bond::Up, "/"),
            (Bond::Down, "\\"),
        ];

        for (bond, expected) in cases {
            assert_eq!(expected, bond.to_string());
        }
    }

    #[test]
    fn directional_bonds_collapse_to_single() {
        assert_eq!(Bond::Up.without_direction(), Bond::Single);
        assert_eq!(Bond::Down.without_direction(), Bond::Single);
        assert_eq!(Bond::Double.without_direction(), Bond::Double);
        assert_eq!(Bond::Aromatic.without_direction(), Bond::Aromatic);
    }
}
