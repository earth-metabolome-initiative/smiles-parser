use alloc::{
    string::{String, ToString},
    vec::Vec,
};
use core::str::FromStr;

use elements_rs::{Element, Isotope};
use molecular_formulas::{ChemicalFormula, errors::ParserError};
use thiserror::Error;

use super::Smiles;

/// Error raised while converting a [`Smiles`] graph into a molecular formula.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum MolecularFormulaConversionError {
    /// The SMILES graph does not contain any atoms.
    #[error("cannot convert an empty SMILES graph to a molecular formula")]
    EmptySmiles,
    /// A wildcard atom cannot be represented by an exact molecular formula.
    #[error("cannot convert wildcard atom at index {atom_id} to a molecular formula")]
    WildcardAtom {
        /// Index of the wildcard atom in the SMILES graph.
        atom_id: usize,
    },
    /// The target molecular formula type cannot represent this isotope.
    #[error("cannot convert unsupported isotope {mass_number}{element} at atom index {atom_id}")]
    UnsupportedIsotope {
        /// Index of the isotope atom in the SMILES graph.
        atom_id: usize,
        /// Element carrying the unsupported isotope mass number.
        element: Element,
        /// Parsed isotope mass number.
        mass_number: u16,
    },
    /// An atom or hydrogen count exceeded the supported conversion range.
    #[error("SMILES molecular formula count overflowed")]
    CountOverflow,
    /// A component charge exceeded the supported conversion range.
    #[error("SMILES molecular formula charge overflowed")]
    ChargeOverflow,
    /// The generated formula was rejected by `molecular-formulas`.
    #[error("generated molecular formula `{formula}` failed to parse: {source}")]
    Parser {
        /// Generated formula string.
        formula: String,
        /// Parser error returned by `molecular-formulas`.
        source: ParserError,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FormulaSpecies {
    Element(Element),
    Isotope { element: Element, mass_number: u16 },
}

impl FormulaSpecies {
    #[inline]
    fn element(self) -> Element {
        match self {
            Self::Element(element) | Self::Isotope { element, .. } => element,
        }
    }

    #[inline]
    fn isotope_mass_number(self) -> Option<u16> {
        match self {
            Self::Element(_) => None,
            Self::Isotope { mass_number, .. } => Some(mass_number),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct FormulaSpeciesCount {
    species: FormulaSpecies,
    count: u32,
}

#[derive(Debug, Clone, Default)]
struct ComponentFormula {
    species_counts: Vec<FormulaSpeciesCount>,
    charge: i32,
}

impl ComponentFormula {
    fn add_species(
        &mut self,
        species: FormulaSpecies,
        count: u32,
    ) -> Result<(), MolecularFormulaConversionError> {
        if count == 0 {
            return Ok(());
        }

        if let Some(entry) = self.species_counts.iter_mut().find(|entry| entry.species == species) {
            entry.count = entry
                .count
                .checked_add(count)
                .ok_or(MolecularFormulaConversionError::CountOverflow)?;
        } else {
            self.species_counts.push(FormulaSpeciesCount { species, count });
        }
        Ok(())
    }

    fn write_formula(&mut self, target: &mut String) {
        self.sort_species();
        for entry in &self.species_counts {
            write_species(target, entry.species);
            if entry.count != 1 {
                target.push_str(&entry.count.to_string());
            }
        }
        write_charge(target, self.charge);
    }

    fn sort_species(&mut self) {
        let has_carbon =
            self.species_counts.iter().any(|entry| entry.species.element() == Element::C);
        self.species_counts.sort_by(|left, right| {
            formula_species_order(left.species, has_carbon)
                .cmp(&formula_species_order(right.species, has_carbon))
        });
    }
}

impl TryFrom<&Smiles> for ChemicalFormula {
    type Error = MolecularFormulaConversionError;

    fn try_from(smiles: &Smiles) -> Result<Self, Self::Error> {
        let formula = smiles_formula_string(smiles)?;
        Self::from_str(&formula)
            .map_err(|source| MolecularFormulaConversionError::Parser { formula, source })
    }
}

impl TryFrom<Smiles> for ChemicalFormula {
    type Error = MolecularFormulaConversionError;

    fn try_from(smiles: Smiles) -> Result<Self, Self::Error> {
        Self::try_from(&smiles)
    }
}

fn smiles_formula_string(smiles: &Smiles) -> Result<String, MolecularFormulaConversionError> {
    if smiles.nodes().is_empty() {
        return Err(MolecularFormulaConversionError::EmptySmiles);
    }

    let components = smiles.connected_components();
    let mut component_formulas =
        vec![ComponentFormula::default(); components.number_of_components()];
    let component_ids = components.component_identifiers().collect::<Vec<_>>();

    for (atom_id, atom) in smiles.nodes().iter().enumerate() {
        let component = &mut component_formulas[component_ids[atom_id]];
        let element =
            atom.element().ok_or(MolecularFormulaConversionError::WildcardAtom { atom_id })?;
        if let Some(mass_number) = atom.isotope_mass_number()
            && Isotope::try_from((element, mass_number)).is_err()
        {
            return Err(MolecularFormulaConversionError::UnsupportedIsotope {
                atom_id,
                element,
                mass_number,
            });
        }
        let species = atom.isotope_mass_number().map_or(FormulaSpecies::Element(element), |mass| {
            FormulaSpecies::Isotope { element, mass_number: mass }
        });
        component.add_species(species, 1)?;
        let hydrogen_count = u32::from(atom.hydrogen_count())
            .checked_add(u32::from(smiles.implicit_hydrogen_count(atom_id)))
            .ok_or(MolecularFormulaConversionError::CountOverflow)?;
        component.add_species(FormulaSpecies::Element(Element::H), hydrogen_count)?;
        component.charge = component
            .charge
            .checked_add(i32::from(atom.charge_value()))
            .ok_or(MolecularFormulaConversionError::ChargeOverflow)?;
    }

    let mut formula = String::new();
    for (index, component) in component_formulas.iter_mut().enumerate() {
        if index != 0 {
            formula.push('.');
        }
        component.write_formula(&mut formula);
    }
    Ok(formula)
}

fn formula_species_order(
    species: FormulaSpecies,
    has_carbon: bool,
) -> (u8, &'static str, u8, bool, Option<u16>) {
    let element = species.element();
    let isotope_mass_number = species.isotope_mass_number();
    let is_isotope = isotope_mass_number.is_some();
    let atomic_number = u8::from(element);

    if has_carbon {
        match (element, is_isotope) {
            (Element::C, _) => (0, "", atomic_number, is_isotope, isotope_mass_number),
            (Element::H, _) => (1, "", atomic_number, is_isotope, isotope_mass_number),
            (_, false) => (2, element.symbol(), 0, false, None),
            (_, true) => (3, "", atomic_number, true, isotope_mass_number),
        }
    } else {
        match (element, is_isotope) {
            (Element::H, _) => (0, "", atomic_number, is_isotope, isotope_mass_number),
            (_, false) => (1, element.symbol(), 0, false, None),
            (_, true) => (2, "", atomic_number, true, isotope_mass_number),
        }
    }
}

fn write_species(target: &mut String, species: FormulaSpecies) {
    match species {
        FormulaSpecies::Element(element) => target.push_str(element.symbol()),
        FormulaSpecies::Isotope { element, mass_number } => {
            target.push('[');
            target.push_str(&mass_number.to_string());
            target.push_str(element.symbol());
            target.push(']');
        }
    }
}

fn write_charge(target: &mut String, charge: i32) {
    match charge.cmp(&0) {
        core::cmp::Ordering::Less => {
            target.push('-');
            if charge != -1 {
                target.push_str(&charge.unsigned_abs().to_string());
            }
        }
        core::cmp::Ordering::Equal => {}
        core::cmp::Ordering::Greater => {
            target.push('+');
            if charge != 1 {
                target.push_str(&charge.to_string());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;

    use molecular_formulas::{ChargedMolecularFormula, MolecularFormula};

    use super::*;

    #[test]
    fn formula_string_counts_implicit_hydrogens() {
        let smiles: Smiles = "c1ccccc1".parse().unwrap();

        assert_eq!(smiles_formula_string(&smiles).unwrap(), "C6H6");
    }

    #[test]
    fn formula_conversion_preserves_isotopes_explicit_hydrogens_and_charge() {
        let smiles: Smiles = "[13CH3][NH3+]".parse().unwrap();
        let formula = ChemicalFormula::try_from(&smiles).unwrap();

        assert_eq!(formula.to_string(), "[¹³C]H₆N⁺");
        assert_eq!(formula.count_of_element::<u32>(Element::C), Some(1));
        assert_eq!(formula.count_of_element::<u32>(Element::H), Some(6));
        assert_eq!(formula.count_of_element::<u32>(Element::N), Some(1));
        assert!((formula.charge() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn formula_conversion_preserves_disconnected_components() {
        let smiles: Smiles = "[Na+].[Cl-]".parse().unwrap();
        let formula = ChemicalFormula::try_from(smiles).unwrap();

        assert_eq!(formula.to_string(), "Na⁺.Cl⁻");
        assert!(formula.charge().abs() < f64::EPSILON);
    }

    #[test]
    fn formula_conversion_matches_rdkit_non_carbon_isotope_ordering() {
        let smiles: Smiles = "C([18OH])([131I])Cl".parse().unwrap();
        let formula = ChemicalFormula::try_from(smiles).unwrap();
        let rdkit_formula = ChemicalFormula::try_from("CH2Cl[18O][131I]").unwrap();

        assert_eq!(formula, rdkit_formula);
    }

    #[test]
    fn formula_conversion_rejects_empty_smiles() {
        let error = ChemicalFormula::try_from(Smiles::new()).unwrap_err();

        assert_eq!(error, MolecularFormulaConversionError::EmptySmiles);
    }

    #[test]
    fn formula_conversion_rejects_wildcards() {
        let smiles: Smiles = "*".parse().unwrap();
        let error = ChemicalFormula::try_from(smiles).unwrap_err();

        assert_eq!(error, MolecularFormulaConversionError::WildcardAtom { atom_id: 0 });
    }

    #[test]
    fn formula_conversion_rejects_unsupported_isotopes() {
        let smiles: Smiles = "[999C]".parse().unwrap();
        let error = ChemicalFormula::try_from(smiles).unwrap_err();

        assert_eq!(
            error,
            MolecularFormulaConversionError::UnsupportedIsotope {
                atom_id: 0,
                element: Element::C,
                mass_number: 999,
            }
        );
    }
}
