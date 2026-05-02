use alloc::vec::Vec;

use geometric_traits::traits::SparseValuedMatrix2DRef;

use super::Smiles;
use crate::{
    atom::{AtomSyntax, atom_symbol::AtomSymbol, bracketed::chirality::Chirality},
    bond::Bond,
};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub(crate) struct BondKindHistogram {
    counts: [usize; 5],
}

impl BondKindHistogram {
    #[inline]
    pub(crate) fn record(&mut self, bond: Bond) {
        self.counts[bond_kind_index(bond)] += 1;
    }

    #[inline]
    #[must_use]
    pub(crate) fn count(self, bond: Bond) -> usize {
        self.counts[bond_kind_index(bond)]
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub(crate) struct AtomInvariant {
    pub(crate) syntax: AtomSyntax,
    pub(crate) symbol: AtomSymbol,
    pub(crate) isotope_mass_number: Option<u16>,
    pub(crate) aromatic: bool,
    pub(crate) hydrogens: u8,
    pub(crate) charge: i8,
    pub(crate) class: u16,
    pub(crate) chirality: Option<Chirality>,
    pub(crate) degree: usize,
    pub(crate) bond_kind_histogram: BondKindHistogram,
}

impl<AtomPolicy: crate::smiles::SmilesAtomPolicy> Smiles<AtomPolicy> {
    #[inline]
    #[must_use]
    pub(crate) fn atom_invariant(&self, node_id: usize) -> Option<AtomInvariant> {
        let atom = *self.node_by_id(node_id)?;
        let mut bond_kind_histogram = BondKindHistogram::default();
        let mut degree = 0;
        for entry in self.bond_matrix.sparse_row_values_ref(node_id) {
            degree += 1;
            bond_kind_histogram.record(entry.bond());
        }

        Some(AtomInvariant {
            syntax: atom.syntax(),
            symbol: atom.symbol(),
            isotope_mass_number: atom.isotope_mass_number(),
            aromatic: atom.aromatic(),
            hydrogens: atom.hydrogen_count(),
            charge: atom.charge_value(),
            class: atom.class(),
            chirality: atom.chirality(),
            degree,
            bond_kind_histogram,
        })
    }

    #[inline]
    #[must_use]
    pub(crate) fn atom_invariants(&self) -> Vec<AtomInvariant> {
        (0..self.atom_nodes.len()).map(|node_id| self.atom_invariant(node_id).unwrap()).collect()
    }
}

#[inline]
pub(crate) const fn bond_kind_index(bond: Bond) -> usize {
    bond_kind_code(bond) as usize
}

#[inline]
pub(crate) const fn bond_kind_code(bond: Bond) -> u8 {
    match bond {
        Bond::Single | Bond::Up | Bond::Down => 0,
        Bond::Double => 1,
        Bond::Triple => 2,
        Bond::Quadruple => 3,
        Bond::Aromatic => 4,
    }
}

#[inline]
pub(crate) const fn planning_chirality_key(chirality: Option<Chirality>) -> (u8, u8) {
    match chirality {
        None => (0, 0),
        Some(Chirality::At | Chirality::AtAt) => (1, 0),
        Some(Chirality::TH(_)) => (2, 0),
        Some(Chirality::AL(_)) => (3, 0),
        Some(Chirality::SP(_)) => (4, 0),
        Some(Chirality::TB(_)) => (5, 0),
        Some(Chirality::OH(_)) => (6, 0),
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use elements_rs::Element;

    use super::{
        AtomInvariant, BondKindHistogram, Smiles, bond_kind_code, bond_kind_index,
        planning_chirality_key,
    };
    use crate::{
        atom::{
            Atom, AtomSyntax,
            atom_symbol::AtomSymbol,
            bracketed::{charge::Charge, chirality::Chirality},
        },
        bond::{
            Bond,
            bond_edge::{BondEdge, bond_edge},
        },
        smiles::BondMatrixBuilder,
    };

    fn smiles_from_edges(atom_nodes: Vec<Atom>, bond_edges: &[BondEdge]) -> Smiles {
        let mut builder = BondMatrixBuilder::with_capacity(bond_edges.len());
        for edge in bond_edges {
            builder.push_edge(edge.0, edge.1, edge.2, edge.3).unwrap();
        }
        let number_of_nodes = atom_nodes.len();
        Smiles::from_bond_matrix_parts(atom_nodes, builder.finish(number_of_nodes))
    }

    #[test]
    fn atom_invariants_of_empty_graph_are_empty() {
        assert!(
            Smiles::<crate::smiles::ConcreteAtoms>::new_for_policy().atom_invariants().is_empty()
        );
    }

    #[test]
    fn atom_invariant_captures_local_atom_fields() {
        let atom = Atom::builder()
            .with_symbol(AtomSymbol::Element(Element::N))
            .with_isotope(15)
            .with_aromatic(true)
            .with_hydrogens(2)
            .with_charge(Charge::try_new(2).unwrap())
            .with_class(17)
            .with_chirality(Chirality::TH(2))
            .build();
        let smiles = smiles_from_edges(vec![atom], &[]);

        let invariant = smiles.atom_invariant(0).unwrap();
        assert_eq!(
            invariant,
            AtomInvariant {
                syntax: AtomSyntax::Bracket,
                symbol: AtomSymbol::Element(Element::N),
                isotope_mass_number: Some(15),
                aromatic: true,
                hydrogens: 2,
                charge: 2,
                class: 17,
                chirality: Some(Chirality::TH(2)),
                degree: 0,
                bond_kind_histogram: BondKindHistogram::default(),
            }
        );
    }

    #[test]
    fn atom_invariants_distinguish_organic_subset_bracket_and_wildcard_atoms() {
        let organic = Atom::new_organic_subset(AtomSymbol::Element(Element::C), false);
        let bracket = Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).build();
        let wildcard = Atom::builder().with_symbol(AtomSymbol::WildCard).build();
        let smiles = smiles_from_edges(vec![organic, bracket, wildcard], &[]);

        let invariants = smiles.atom_invariants();
        assert_ne!(invariants[0], invariants[1]);
        assert_ne!(invariants[1], invariants[2]);
        assert_eq!(invariants[0].syntax, AtomSyntax::OrganicSubset);
        assert_eq!(invariants[1].syntax, AtomSyntax::Bracket);
        assert_eq!(invariants[2].symbol, AtomSymbol::WildCard);
    }

    #[test]
    fn atom_invariant_captures_degree_and_incident_bond_kinds() {
        let smiles = smiles_from_edges(
            vec![
                Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::O), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::N), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::S), true),
            ],
            &[
                bond_edge(0, 1, Bond::Single, None),
                bond_edge(0, 2, Bond::Double, None),
                bond_edge(0, 3, Bond::Aromatic, None),
            ],
        );

        let invariant = smiles.atom_invariant(0).unwrap();
        assert_eq!(invariant.degree, 3);
        assert_eq!(invariant.bond_kind_histogram.count(Bond::Single), 1);
        assert_eq!(invariant.bond_kind_histogram.count(Bond::Double), 1);
        assert_eq!(invariant.bond_kind_histogram.count(Bond::Aromatic), 1);
    }

    #[test]
    fn atom_invariants_distinguish_equal_degree_but_different_bond_histograms() {
        let single_single = smiles_from_edges(
            vec![
                Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::O), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::N), false),
            ],
            &[bond_edge(0, 1, Bond::Single, None), bond_edge(0, 2, Bond::Single, None)],
        );
        let single_double = smiles_from_edges(
            vec![
                Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::O), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::N), false),
            ],
            &[bond_edge(0, 1, Bond::Single, None), bond_edge(0, 2, Bond::Double, None)],
        );

        assert_ne!(single_single.atom_invariant(0), single_double.atom_invariant(0));
    }

    #[test]
    fn atom_invariant_treats_directional_single_bonds_as_single_family() {
        let smiles = smiles_from_edges(
            vec![
                Atom::new_organic_subset(AtomSymbol::Element(Element::C), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::O), false),
                Atom::new_organic_subset(AtomSymbol::Element(Element::N), false),
            ],
            &[bond_edge(0, 1, Bond::Up, None), bond_edge(0, 2, Bond::Down, None)],
        );

        let invariant = smiles.atom_invariant(0).unwrap();
        assert_eq!(invariant.bond_kind_histogram.count(Bond::Single), 2);
        assert_eq!(invariant.bond_kind_histogram.count(Bond::Up), 2);
        assert_eq!(invariant.bond_kind_histogram.count(Bond::Down), 2);
    }

    #[test]
    fn bond_kind_helpers_group_directional_single_bonds_and_separate_others() {
        assert_eq!(bond_kind_code(Bond::Single), 0);
        assert_eq!(bond_kind_code(Bond::Up), 0);
        assert_eq!(bond_kind_code(Bond::Down), 0);
        assert_eq!(bond_kind_code(Bond::Double), 1);
        assert_eq!(bond_kind_code(Bond::Triple), 2);
        assert_eq!(bond_kind_code(Bond::Quadruple), 3);
        assert_eq!(bond_kind_code(Bond::Aromatic), 4);

        assert_eq!(bond_kind_index(Bond::Single), bond_kind_index(Bond::Up));
        assert_eq!(bond_kind_index(Bond::Up), bond_kind_index(Bond::Down));
        assert_ne!(bond_kind_index(Bond::Single), bond_kind_index(Bond::Double));
        assert_ne!(bond_kind_index(Bond::Triple), bond_kind_index(Bond::Aromatic));
    }

    #[test]
    fn planning_chirality_key_groups_surface_variants_by_family() {
        assert_eq!(planning_chirality_key(None), (0, 0));
        assert_eq!(planning_chirality_key(Some(Chirality::At)), (1, 0));
        assert_eq!(planning_chirality_key(Some(Chirality::AtAt)), (1, 0));
        assert_eq!(planning_chirality_key(Some(Chirality::TH(1))), (2, 0));
        assert_eq!(planning_chirality_key(Some(Chirality::AL(2))), (3, 0));
        assert_eq!(planning_chirality_key(Some(Chirality::SP(3))), (4, 0));
        assert_eq!(planning_chirality_key(Some(Chirality::TB(20))), (5, 0));
        assert_eq!(planning_chirality_key(Some(Chirality::OH(30))), (6, 0));
    }
}
