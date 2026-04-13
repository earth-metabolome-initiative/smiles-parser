use alloc::vec::Vec;

use elements_rs::Element;
use geometric_traits::{
    prelude::canonical_label_labeled_simple_graph,
    traits::{SparseMatrix2D, SparseValuedMatrix2DRef, SparseValuedMatrixRef},
};

use super::{
    BondMatrix, BondMatrixBuilder, Smiles, StereoNeighbor,
    double_bond_stereo::DoubleBondStereoConfig,
    implicit_hydrogens::implicit_hydrogens_if_written_unbracketed,
    stereo::normalized_tetrahedral_chirality,
};
use crate::{
    atom::{Atom, AtomSyntax, atom_symbol::AtomSymbol, bracketed::chirality::Chirality},
    bond::Bond,
};

#[cfg(any(test, feature = "fuzzing"))]
mod support;
#[cfg(test)]
mod tests;

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct CanonicalAtomLabel {
    syntax: u8,
    symbol: AtomSymbol,
    isotope_mass_number: Option<u16>,
    aromatic: bool,
    hydrogens: u8,
    charge: i8,
    class: u16,
    chirality_kind: u8,
    chirality_value: u8,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct CanonicalBondLabel(u8);

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct CanonicalStereoNeighborKey(u8, usize);

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct CanonicalizationStateKey {
    atom_labels: Vec<CanonicalAtomLabel>,
    bond_edges: Vec<(usize, usize, CanonicalBondLabel)>,
    parsed_stereo_neighbors: Vec<Vec<CanonicalStereoNeighborKey>>,
    implicit_hydrogen_cache: Option<Vec<u8>>,
}

#[derive(Debug, Clone)]
struct AtomBasedDoubleBondNormalization {
    override_rows: Vec<Vec<(usize, Bond)>>,
    semantic_endpoints: Vec<bool>,
    clear_chirality: Vec<bool>,
}

#[derive(Debug, Clone)]
struct NonSemanticDirectionalNormalization {
    override_rows: Vec<Vec<(usize, Bond)>>,
}

#[derive(Debug, Clone, Copy)]
struct AtomBasedDoubleBondSide {
    endpoint: usize,
    reference_atom: usize,
    reference_bond_is_up: bool,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum StereoSubstituentIdentityKey {
    ExplicitHydrogen,
    Atom(AtomBasedSubstituentPriorityKey),
}

/// Canonical labeling result for a [`Smiles`] graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesCanonicalLabeling {
    order: Vec<usize>,
    new_index_of_old_node: Vec<usize>,
}

impl SmilesCanonicalLabeling {
    fn new(order: Vec<usize>) -> Self {
        let mut new_index_of_old_node = vec![usize::MAX; order.len()];
        for (new_index, old_node) in order.iter().copied().enumerate() {
            new_index_of_old_node[old_node] = new_index;
        }
        Self { order, new_index_of_old_node }
    }

    /// Returns original node ids in canonical order.
    #[inline]
    #[must_use]
    pub fn order(&self) -> &[usize] {
        &self.order
    }

    /// Returns the new canonical index for each original node id.
    #[inline]
    #[must_use]
    pub fn new_index_of_old_node(&self) -> &[usize] {
        &self.new_index_of_old_node
    }
}

impl Smiles {
    fn exact_canonical_labeling(&self) -> SmilesCanonicalLabeling {
        let result = canonical_label_labeled_simple_graph(
            self,
            |node_id| canonical_atom_label(self.nodes()[node_id]),
            |node_a, node_b| {
                canonical_bond_label(
                    self.edge_for_node_pair((node_a, node_b))
                        .unwrap_or_else(|| unreachable!("canonizer only queries existing edges"))
                        .bond(),
                )
            },
        );
        SmilesCanonicalLabeling::new(result.order)
    }

    fn stereo_neutral_canonical_labeling(&self) -> SmilesCanonicalLabeling {
        let result = canonical_label_labeled_simple_graph(
            self,
            |node_id| stereo_neutral_canonical_atom_label(self.nodes()[node_id]),
            |node_a, node_b| {
                canonical_bond_label(
                    self.edge_for_node_pair((node_a, node_b))
                        .unwrap_or_else(|| unreachable!("canonizer only queries existing edges"))
                        .bond()
                        .without_direction(),
                )
            },
        );
        SmilesCanonicalLabeling::new(result.order)
    }

    fn exact_canonicalize(&self) -> Self {
        let labeling = self.exact_canonical_labeling();
        let order = labeling.order();
        let new_index_of_old_node = labeling.new_index_of_old_node();

        let atom_nodes: Vec<_> =
            order.iter().copied().map(|old_node| self.atom_nodes[old_node]).collect();

        let mut canonical_edges = Vec::with_capacity(self.number_of_bonds());
        for ((row, column), entry) in self.bond_matrix().sparse_entries() {
            if row >= column {
                continue;
            }
            let new_row = new_index_of_old_node[row];
            let new_column = new_index_of_old_node[column];
            let (new_row, new_column) = Smiles::edge_key(new_row, new_column);
            canonical_edges.push((new_row, new_column, entry.bond()));
        }
        canonical_edges.sort_unstable_by_key(|(row, column, _bond)| (*row, *column));
        let bond_matrix = BondMatrix::from_sorted_upper_triangular_entries(
            atom_nodes.len(),
            canonical_edges.into_iter().enumerate().map(
                |(canonical_order, (row, column, bond))| {
                    (row, column, super::BondEntry::new(bond, None, canonical_order))
                },
            ),
        )
        .unwrap_or_else(|_| {
            unreachable!("canonical relabeling preserves a simple upper-triangular graph")
        });

        let parsed_stereo_neighbors = order
            .iter()
            .copied()
            .map(|old_node| {
                remap_parsed_stereo_neighbors_row(self, old_node, new_index_of_old_node)
            })
            .collect();

        let implicit_hydrogen_cache = self
            .implicit_hydrogen_cache
            .as_ref()
            .map(|cache| order.iter().copied().map(|old_node| cache[old_node]).collect());

        Self::from_bond_matrix_parts_with_sidecars(
            atom_nodes,
            bond_matrix,
            parsed_stereo_neighbors,
            implicit_hydrogen_cache,
            None,
        )
    }

    fn stereo_normal_form(&self) -> Self {
        if self.nodes().is_empty() {
            return self.clone();
        }

        let stereo_neutral_labeling = self.stereo_neutral_canonical_labeling();
        let new_index_of_old_node = stereo_neutral_labeling.new_index_of_old_node();
        let refined_classes = self.stereo_neutral_refined_classes();
        let rooted_classes = self.stereo_neutral_rooted_classes(&refined_classes);
        let directional_overrides = self.projected_directional_bond_overrides_with_classes(
            new_index_of_old_node,
            &rooted_classes,
            &refined_classes,
        );
        let atom_based_double_bond_normalization = self.atom_based_double_bond_normalization(
            new_index_of_old_node,
            &rooted_classes,
            &refined_classes,
        );
        let non_semantic_directional_normalization =
            self.non_semantic_directional_normalization(&rooted_classes, &refined_classes);

        let normalized_atom_rows = self.stereo_normalized_atom_rows(
            new_index_of_old_node,
            &rooted_classes,
            &refined_classes,
            &atom_based_double_bond_normalization,
        );

        let atom_nodes = normalized_atom_rows.iter().map(|(atom, _)| *atom).collect::<Vec<_>>();
        let parsed_stereo_neighbors = normalized_atom_rows
            .into_iter()
            .map(|(_atom, neighbors)| neighbors)
            .collect::<Vec<_>>();

        let mut builder = BondMatrixBuilder::with_capacity(self.number_of_bonds());
        for ((row, column), entry) in self.bond_matrix().sparse_entries() {
            if row >= column {
                continue;
            }

            let mut bond = entry.bond();
            if let Some(override_bond) = atom_based_override_bond(
                &non_semantic_directional_normalization.override_rows,
                row,
                column,
            ) {
                bond = override_bond;
            }
            if matches!(bond, Bond::Up | Bond::Down)
                && (directional_overrides.has_semantic_endpoint(row)
                    || directional_overrides.has_semantic_endpoint(column)
                    || atom_based_double_bond_normalization.semantic_endpoints[row]
                    || atom_based_double_bond_normalization.semantic_endpoints[column])
            {
                bond = Bond::Single;
            }
            if let Some(override_bond) = directional_overrides.get(row, column) {
                bond = override_bond;
            }
            if let Some(override_bond) = atom_based_override_bond(
                &atom_based_double_bond_normalization.override_rows,
                row,
                column,
            ) {
                bond = override_bond;
            }

            builder
                .push_edge(row, column, bond, None)
                .unwrap_or_else(|_| unreachable!("stereo normalization preserves a simple graph"));
        }

        Self::from_bond_matrix_parts_with_sidecars(
            atom_nodes,
            builder.finish(self.nodes().len()),
            parsed_stereo_neighbors,
            None,
            None,
        )
    }

    fn stereo_normalized_atom_rows(
        &self,
        new_index_of_old_node: &[usize],
        rooted_classes: &[usize],
        refined_classes: &[usize],
        atom_based_double_bond_normalization: &AtomBasedDoubleBondNormalization,
    ) -> Vec<(Atom, Vec<StereoNeighbor>)> {
        self.atom_nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(node_id, atom)| {
                if atom_based_double_bond_normalization.clear_chirality[node_id] {
                    let cleared = atom_with_chirality(atom, None);
                    return (
                        maybe_collapse_atom_to_organic_subset(self, node_id, cleared),
                        Vec::new(),
                    );
                }
                let parsed_neighbors = self.parsed_stereo_neighbors(node_id);
                let normalized_neighbors = canonical_stereo_neighbors_row(
                    self,
                    node_id,
                    atom.chirality(),
                    &parsed_neighbors,
                    new_index_of_old_node,
                );
                let normalized_chirality = normalized_stereo_chirality(
                    self,
                    node_id,
                    atom.chirality(),
                    &parsed_neighbors,
                    &normalized_neighbors,
                    rooted_classes,
                    refined_classes,
                );
                let normalized_atom = atom_with_chirality(atom, normalized_chirality);
                (
                    if atom.chirality().is_some() && normalized_chirality.is_none() {
                        maybe_collapse_atom_to_organic_subset(self, node_id, normalized_atom)
                    } else {
                        normalized_atom
                    },
                    normalized_chirality.map_or_else(Vec::new, |_| normalized_neighbors),
                )
            })
            .collect::<Vec<_>>()
    }

    fn canonicalization_normal_form(&self) -> Self {
        let input = self.canonicalization_spelling_normal_form();
        if input.nodes().is_empty() {
            return input;
        }
        let Ok(perception) = input.perceive_aromaticity() else {
            return input;
        };
        perception.into_aromaticized().canonicalization_spelling_normal_form()
    }

    fn canonicalization_step(&self) -> Self {
        let canonicalized = self
            .canonicalization_normal_form()
            .collapse_removable_explicit_hydrogens()
            .stereo_normal_form()
            .exact_canonicalize()
            .canonicalization_spelling_normal_form();
        let has_aromatic_bonds = canonicalized
            .bond_matrix()
            .sparse_entries()
            .any(|((_row, _column), entry)| entry.bond() == Bond::Aromatic);
        if !has_aromatic_bonds {
            return canonicalized;
        }

        canonicalized
            .kekulize_standalone()
            .ok()
            .map(|kekulized| {
                kekulized
                    .stereo_normal_form()
                    .exact_canonicalize()
                    .canonicalization_spelling_normal_form()
            })
            .unwrap_or(canonicalized)
    }

    fn canonicalization_spelling_normal_form(&self) -> Self {
        let atom_nodes = self
            .atom_nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(node_id, atom)| canonicalization_atom_spelling_normal_form(self, node_id, atom))
            .collect();
        Self::from_bond_matrix_parts_with_sidecars(
            atom_nodes,
            self.bond_matrix.clone(),
            self.parsed_stereo_neighbors.clone(),
            self.implicit_hydrogen_cache.clone(),
            None,
        )
    }

    /// Returns the canonical labeling of the current graph.
    #[must_use]
    pub fn canonical_labeling(&self) -> SmilesCanonicalLabeling {
        self.canonicalization_normal_form()
            .collapse_removable_explicit_hydrogens()
            .stereo_normal_form()
            .exact_canonical_labeling()
    }

    /// Returns whether the current graph is already in canonical form.
    #[must_use]
    pub fn is_canonical(&self) -> bool {
        let canonicalized = self.canonicalize();
        self.atom_nodes == canonicalized.atom_nodes
            && self.bond_matrix == canonicalized.bond_matrix
            && self.parsed_stereo_neighbors == canonicalized.parsed_stereo_neighbors
            && self.implicit_hydrogen_cache == canonicalized.implicit_hydrogen_cache
            && self.kekulization_source == canonicalized.kekulization_source
    }

    /// Returns the graph rewritten into canonical node order.
    #[must_use]
    pub fn canonicalize(&self) -> Self {
        self.canonicalize_orbit_min()
    }

    fn canonicalize_orbit_min(&self) -> Self {
        let mut states: Vec<Self> = Vec::new();
        let mut keys: Vec<CanonicalizationStateKey> = Vec::new();
        let mut current = self.canonicalization_step();

        loop {
            let key = canonicalization_state_key(&current);
            if let Some(cycle_start) = keys.iter().position(|existing| *existing == key) {
                let best_relative_index = keys[cycle_start..]
                    .iter()
                    .enumerate()
                    .min_by(|left, right| left.1.cmp(right.1))
                    .map_or_else(|| unreachable!("cycle slice is non-empty"), |(index, _)| index);
                return states[cycle_start + best_relative_index].clone();
            }
            keys.push(key);
            states.push(current.clone());
            current = current.canonicalization_step();
        }
    }

    #[allow(clippy::too_many_lines)]
    fn atom_based_double_bond_normalization(
        &self,
        preorder_indices: &[usize],
        rooted_classes: &[usize],
        refined_classes: &[usize],
    ) -> AtomBasedDoubleBondNormalization {
        let mut records = Vec::new();
        let mut clear_chirality = vec![false; self.nodes().len()];

        for ((row, column), entry) in self.bond_matrix().sparse_entries() {
            if row >= column || entry.bond() != Bond::Double {
                continue;
            }
            if !self.atom_based_double_bond_supports_semantic_stereo(row, column) {
                continue;
            }
            let Some((side_a, side_b)) =
                self.atom_based_double_bond_sides(row, column, rooted_classes, refined_classes)
            else {
                continue;
            };
            let config = if side_a.reference_bond_is_up == side_b.reference_bond_is_up {
                DoubleBondStereoConfig::E
            } else {
                DoubleBondStereoConfig::Z
            };
            records.push((row, column, side_a, side_b, config));
            clear_chirality[row] = true;
            clear_chirality[column] = true;
        }

        let mut semantic_endpoints = vec![false; self.nodes().len()];
        let mut rows = vec![Vec::new(); self.nodes().len()];
        if records.is_empty() {
            return AtomBasedDoubleBondNormalization {
                override_rows: rows,
                semantic_endpoints,
                clear_chirality,
            };
        }

        let mut edge_keys = records
            .iter()
            .flat_map(|record| {
                [
                    Smiles::edge_key(record.2.endpoint, record.2.reference_atom),
                    Smiles::edge_key(record.3.endpoint, record.3.reference_atom),
                ]
            })
            .collect::<Vec<_>>();
        edge_keys.sort_unstable();
        edge_keys.dedup();
        let mut adjacency: Vec<Vec<(usize, bool)>> = vec![Vec::new(); edge_keys.len()];

        for (_row, _column, side_a, side_b, config) in &records {
            semantic_endpoints[side_a.endpoint] = true;
            semantic_endpoints[side_b.endpoint] = true;
            let left_edge = edge_keys
                .binary_search(&Smiles::edge_key(side_a.endpoint, side_a.reference_atom))
                .unwrap_or_else(|_| unreachable!());
            let right_edge = edge_keys
                .binary_search(&Smiles::edge_key(side_b.endpoint, side_b.reference_atom))
                .unwrap_or_else(|_| unreachable!());
            let same_parity = matches!(config, DoubleBondStereoConfig::E);
            adjacency[left_edge].push((right_edge, same_parity));
            adjacency[right_edge].push((left_edge, same_parity));
        }

        let mut bond_is_up: Vec<Option<bool>> = vec![None; edge_keys.len()];
        let mut component_seen = vec![false; edge_keys.len()];
        let mut stack = Vec::new();

        for start in 0..edge_keys.len() {
            if component_seen[start] {
                continue;
            }

            let mut component = Vec::new();
            stack.push(start);

            while let Some(current) = stack.pop() {
                if component_seen[current] {
                    continue;
                }
                component_seen[current] = true;
                component.push(current);
                for &(neighbor, _same_parity) in &adjacency[current] {
                    if !component_seen[neighbor] {
                        stack.push(neighbor);
                    }
                }
            }

            let seed = component
                .iter()
                .copied()
                .min_by_key(|&edge_id| {
                    canonical_directional_edge_key(edge_keys[edge_id], preorder_indices)
                })
                .unwrap_or_else(|| unreachable!());
            bond_is_up[seed] = Some(true);
            stack.push(seed);

            while let Some(current) = stack.pop() {
                let current_is_up = bond_is_up[current].unwrap_or_else(|| unreachable!());
                for &(neighbor, same_parity) in &adjacency[current] {
                    let expected = if same_parity { current_is_up } else { !current_is_up };
                    if let Some(existing) = bond_is_up[neighbor] {
                        debug_assert_eq!(existing, expected);
                    } else {
                        bond_is_up[neighbor] = Some(expected);
                        stack.push(neighbor);
                    }
                }
            }
        }

        for (edge_key, is_up) in edge_keys.into_iter().zip(bond_is_up) {
            let bond = if is_up.unwrap_or(true) { Bond::Up } else { Bond::Down };
            rows[edge_key.0].push((edge_key.1, bond));
        }
        for row in &mut rows {
            row.sort_unstable_by_key(|&(neighbor, _)| neighbor);
        }

        AtomBasedDoubleBondNormalization {
            override_rows: rows,
            semantic_endpoints,
            clear_chirality,
        }
    }

    fn non_semantic_directional_normalization(
        &self,
        rooted_classes: &[usize],
        refined_classes: &[usize],
    ) -> NonSemanticDirectionalNormalization {
        let mut rows = vec![Vec::new(); self.nodes().len()];

        for (endpoint_a, endpoint_b) in self
            .bond_matrix()
            .sparse_entries()
            .filter_map(|((row, column), entry)| {
                (row < column && entry.bond() == Bond::Double).then_some((row, column))
            })
            .filter(|&(endpoint_a, endpoint_b)| {
                non_semantic_double_bond_supports_semantic_stereo(self, endpoint_a, endpoint_b)
                    && !non_semantic_double_bond_directional_neighbors(self, endpoint_a, endpoint_b)
                        .is_empty()
                    && !non_semantic_double_bond_directional_neighbors(self, endpoint_b, endpoint_a)
                        .is_empty()
                    && !non_semantic_double_bond_is_in_cycle(self, endpoint_a, endpoint_b)
            })
        {
            let side_a = non_semantic_double_bond_has_unique_reference_substituent(
                self,
                endpoint_a,
                endpoint_b,
                rooted_classes,
                refined_classes,
            );
            let side_b = non_semantic_double_bond_has_unique_reference_substituent(
                self,
                endpoint_b,
                endpoint_a,
                rooted_classes,
                refined_classes,
            );

            if side_a.is_some() && side_b.is_some() {
                continue;
            }

            for directional_neighbor in
                non_semantic_double_bond_directional_neighbors(self, endpoint_a, endpoint_b)
            {
                let (left, right) = Smiles::edge_key(endpoint_a, directional_neighbor);
                rows[left].push((right, Bond::Single));
            }
            for directional_neighbor in
                non_semantic_double_bond_directional_neighbors(self, endpoint_b, endpoint_a)
            {
                let (left, right) = Smiles::edge_key(endpoint_b, directional_neighbor);
                rows[left].push((right, Bond::Single));
            }
        }

        for row in &mut rows {
            row.sort_unstable_by_key(|&(neighbor, _)| neighbor);
            row.dedup_by_key(|entry| entry.0);
        }

        NonSemanticDirectionalNormalization { override_rows: rows }
    }

    fn collapse_removable_explicit_hydrogens(&self) -> Self {
        let node_count = self.nodes().len();
        let mut collapsed_parent_of = vec![None; node_count];
        let mut collapsed_count_for_parent = vec![0_u8; node_count];

        for (node_id, collapsed_parent) in collapsed_parent_of.iter_mut().enumerate() {
            let Some(parent) = self.collapsible_explicit_hydrogen_parent(node_id) else {
                continue;
            };
            *collapsed_parent = Some(parent);
            collapsed_count_for_parent[parent] =
                collapsed_count_for_parent[parent].saturating_add(1);
        }

        if collapsed_parent_of.iter().all(Option::is_none) {
            return self.clone();
        }

        let kept_nodes = (0..node_count)
            .filter(|&node_id| collapsed_parent_of[node_id].is_none())
            .collect::<Vec<_>>();
        let mut new_index_of_old_node = vec![usize::MAX; node_count];
        for (new_index, old_node) in kept_nodes.iter().copied().enumerate() {
            new_index_of_old_node[old_node] = new_index;
        }

        let atom_nodes = kept_nodes
            .iter()
            .copied()
            .map(|old_node| {
                atom_with_hydrogen_count(
                    self.atom_nodes[old_node],
                    self.atom_nodes[old_node]
                        .hydrogen_count()
                        .saturating_add(collapsed_count_for_parent[old_node]),
                )
            })
            .collect::<Vec<_>>();

        let mut builder = BondMatrixBuilder::with_capacity(self.number_of_bonds());
        for ((row, column), entry) in self.bond_matrix().sparse_entries() {
            if row >= column
                || collapsed_parent_of[row].is_some()
                || collapsed_parent_of[column].is_some()
            {
                continue;
            }
            builder
                .push_edge(
                    new_index_of_old_node[row],
                    new_index_of_old_node[column],
                    entry.bond(),
                    None,
                )
                .unwrap_or_else(|_| {
                    unreachable!("collapsing terminal hydrogens preserves simplicity")
                });
        }

        let parsed_stereo_neighbors = kept_nodes
            .iter()
            .copied()
            .map(|old_node| {
                self.parsed_stereo_neighbors(old_node)
                    .into_iter()
                    .filter_map(|neighbor| {
                        match neighbor {
                            StereoNeighbor::ExplicitHydrogen => {
                                Some(StereoNeighbor::ExplicitHydrogen)
                            }
                            StereoNeighbor::Atom(neighbor_id) => {
                                if collapsed_parent_of[neighbor_id] == Some(old_node) {
                                    Some(StereoNeighbor::ExplicitHydrogen)
                                } else if collapsed_parent_of[neighbor_id].is_some() {
                                    None
                                } else {
                                    Some(StereoNeighbor::Atom(new_index_of_old_node[neighbor_id]))
                                }
                            }
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        Self::from_bond_matrix_parts_with_sidecars(
            atom_nodes,
            builder.finish(kept_nodes.len()),
            parsed_stereo_neighbors,
            None,
            None,
        )
    }

    fn collapsible_explicit_hydrogen_parent(&self, node_id: usize) -> Option<usize> {
        let atom = self.nodes()[node_id];
        if atom.element() != Some(Element::H)
            || atom.syntax() != AtomSyntax::Bracket
            || atom.isotope_mass_number().is_some()
            || atom.aromatic()
            || atom.hydrogen_count() != 0
            || atom.charge_value() != 0
            || atom.class() != 0
            || atom.chirality().is_some()
        {
            return None;
        }

        let mut neighbors = self.bond_matrix().sparse_row(node_id);
        let parent = neighbors.next()?;
        if neighbors.next().is_some() {
            return None;
        }

        let bond = self.edge_for_node_pair((node_id, parent))?.bond();
        if !matches!(bond, Bond::Single | Bond::Up | Bond::Down) {
            return None;
        }

        let parent_atom = self.nodes()[parent];
        if parent_atom.element() == Some(Element::H) {
            return None;
        }
        Some(parent)
    }
}

#[cfg(feature = "benchmarking")]
impl Smiles {
    /// Benchmark-only hook for the canonicalization normal-form stage.
    #[doc(hidden)]
    #[must_use]
    pub fn benchmark_canonicalization_normal_form(&self) -> Self {
        self.canonicalization_normal_form()
    }

    /// Benchmark-only hook for removable explicit-hydrogen collapse.
    #[doc(hidden)]
    #[must_use]
    pub fn benchmark_collapse_removable_explicit_hydrogens(&self) -> Self {
        self.collapse_removable_explicit_hydrogens()
    }

    /// Benchmark-only hook for stereo normalization.
    #[doc(hidden)]
    #[must_use]
    pub fn benchmark_stereo_normal_form(&self) -> Self {
        self.stereo_normal_form()
    }

    /// Benchmark-only hook for the exact canonical labeling stage.
    #[doc(hidden)]
    #[must_use]
    pub fn benchmark_exact_canonical_labeling(&self) -> SmilesCanonicalLabeling {
        self.exact_canonical_labeling()
    }

    /// Benchmark-only hook for exact canonical graph rebuilding.
    #[doc(hidden)]
    #[must_use]
    pub fn benchmark_exact_canonicalize(&self) -> Self {
        self.exact_canonicalize()
    }

    /// Benchmark-only hook for one full canonicalization step.
    #[doc(hidden)]
    #[must_use]
    pub fn benchmark_canonicalization_step(&self) -> Self {
        self.canonicalization_step()
    }
}

#[cfg(feature = "fuzzing")]
impl Smiles {
    /// Panics if canonicalization invariants are violated.
    ///
    /// This is intended for fuzzing and other internal validation passes.
    #[doc(hidden)]
    pub fn debug_assert_canonicalization_invariants(&self) {
        support::assert_canonicalization_invariants(self);
    }
}

impl Smiles {
    fn atom_based_double_bond_supports_semantic_stereo(
        &self,
        node_a: usize,
        node_b: usize,
    ) -> bool {
        self.atom_based_non_single_family_bond_count(node_a) == 1
            && self.atom_based_non_single_family_bond_count(node_b) == 1
            && !self.atom_based_double_bond_is_in_cycle(node_a, node_b)
            && self.atom_based_endpoint_supports_semantic_stereo(node_a, node_b)
            && self.atom_based_endpoint_supports_semantic_stereo(node_b, node_a)
    }

    fn atom_based_non_single_family_bond_count(&self, node_id: usize) -> usize {
        self.bond_matrix()
            .sparse_row_values_ref(node_id)
            .filter(|entry| !matches!(entry.bond(), Bond::Single | Bond::Up | Bond::Down))
            .count()
    }

    fn atom_based_double_bond_is_in_cycle(&self, node_a: usize, node_b: usize) -> bool {
        let mut queue = alloc::collections::VecDeque::from([node_a]);
        let mut seen = vec![false; self.nodes().len()];
        seen[node_a] = true;

        while let Some(current) = queue.pop_front() {
            for neighbor in self.bond_matrix().sparse_row(current) {
                if (current == node_a && neighbor == node_b)
                    || (current == node_b && neighbor == node_a)
                {
                    continue;
                }
                if neighbor == node_b {
                    return true;
                }
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        false
    }

    fn atom_based_endpoint_supports_semantic_stereo(
        &self,
        endpoint: usize,
        opposite_endpoint: usize,
    ) -> bool {
        let atom = self.nodes()[endpoint];
        matches!(
            stereo_chirality_normal_form(
                self,
                endpoint,
                atom.chirality(),
                &self.parsed_stereo_neighbors(endpoint),
            ),
            Some(Chirality::At | Chirality::AtAt | Chirality::TH(1 | 2))
        ) && self
            .parsed_stereo_neighbors(endpoint)
            .iter()
            .filter(|neighbor| **neighbor != StereoNeighbor::Atom(opposite_endpoint))
            .count()
            == 2
            && self
                .parsed_stereo_neighbors(endpoint)
                .contains(&StereoNeighbor::Atom(opposite_endpoint))
    }

    fn atom_based_double_bond_sides(
        &self,
        endpoint_a: usize,
        endpoint_b: usize,
        rooted_classes: &[usize],
        refined_classes: &[usize],
    ) -> Option<(AtomBasedDoubleBondSide, AtomBasedDoubleBondSide)> {
        let side_a = self.atom_based_double_bond_side(
            endpoint_a,
            endpoint_b,
            rooted_classes,
            refined_classes,
        )?;
        let side_b = self.atom_based_double_bond_side(
            endpoint_b,
            endpoint_a,
            rooted_classes,
            refined_classes,
        )?;
        Some((side_a, side_b))
    }

    fn atom_based_double_bond_side(
        &self,
        endpoint: usize,
        opposite_endpoint: usize,
        rooted_classes: &[usize],
        refined_classes: &[usize],
    ) -> Option<AtomBasedDoubleBondSide> {
        let parsed_neighbors = self.parsed_stereo_neighbors(endpoint);
        let atom = self.nodes()[endpoint];
        let chirality =
            stereo_chirality_normal_form(self, endpoint, atom.chirality(), &parsed_neighbors)?;
        if !matches!(chirality, Chirality::At | Chirality::AtAt | Chirality::TH(1 | 2)) {
            return None;
        }
        let reference_atom = self.atom_based_reference_substituent(
            endpoint,
            opposite_endpoint,
            rooted_classes,
            refined_classes,
        )?;
        let target_neighbors = atom_based_double_bond_target_neighbors(
            &parsed_neighbors,
            opposite_endpoint,
            reference_atom,
        )?;
        let permutation = atom_based_permutation_from(&parsed_neighbors, &target_neighbors)?;
        let reference_bond_is_up = atom_based_chirality_is_clockwise(chirality)
            ^ atom_based_permutation_is_odd(&permutation);
        Some(AtomBasedDoubleBondSide { endpoint, reference_atom, reference_bond_is_up })
    }

    fn atom_based_reference_substituent(
        &self,
        endpoint: usize,
        opposite_endpoint: usize,
        rooted_classes: &[usize],
        refined_classes: &[usize],
    ) -> Option<usize> {
        let neighbors = self
            .parsed_stereo_neighbors(endpoint)
            .into_iter()
            .filter_map(|neighbor| {
                match neighbor {
                    StereoNeighbor::Atom(node_id) if node_id != opposite_endpoint => Some(node_id),
                    _ => None,
                }
            })
            .collect::<Vec<_>>();
        let (&first, rest) = neighbors.split_first()?;
        let mut best = first;
        let mut best_key = atom_based_substituent_priority_key(
            self,
            endpoint,
            best,
            rooted_classes,
            refined_classes,
        );
        let mut unique_best = true;
        for &candidate in rest {
            let candidate_key = atom_based_substituent_priority_key(
                self,
                endpoint,
                candidate,
                rooted_classes,
                refined_classes,
            );
            match candidate_key.cmp(&best_key) {
                core::cmp::Ordering::Greater => {
                    best = candidate;
                    best_key = candidate_key;
                    unique_best = true;
                }
                core::cmp::Ordering::Equal => unique_best = false,
                core::cmp::Ordering::Less => {}
            }
        }
        unique_best.then_some(best)
    }
}

fn canonical_atom_label(atom: Atom) -> CanonicalAtomLabel {
    let (chirality_kind, chirality_value) = canonical_chirality_key(atom.chirality());

    CanonicalAtomLabel {
        syntax: match atom.syntax() {
            AtomSyntax::OrganicSubset => 0,
            AtomSyntax::Bracket => 1,
        },
        symbol: atom.symbol(),
        isotope_mass_number: atom.isotope_mass_number(),
        aromatic: atom.aromatic(),
        hydrogens: atom.hydrogen_count(),
        charge: atom.charge_value(),
        class: atom.class(),
        chirality_kind,
        chirality_value,
    }
}

fn stereo_neutral_canonical_atom_label(atom: Atom) -> CanonicalAtomLabel {
    let mut label = canonical_atom_label(atom);
    label.chirality_kind = 0;
    label.chirality_value = 0;
    label
}

fn atom_with_chirality(atom: Atom, chirality: Option<Chirality>) -> Atom {
    match atom.syntax() {
        AtomSyntax::OrganicSubset => {
            debug_assert!(
                chirality.is_none(),
                "organic-subset atoms should not carry explicit chirality tags"
            );
            Atom::new_organic_subset(atom.symbol(), atom.aromatic())
        }
        AtomSyntax::Bracket => {
            let mut builder = Atom::builder()
                .with_symbol(atom.symbol())
                .with_aromatic(atom.aromatic())
                .with_hydrogens(atom.hydrogen_count())
                .with_charge(atom.charge())
                .with_class(atom.class());
            if let Some(isotope) = atom.isotope_mass_number() {
                builder = builder.with_isotope(isotope);
            }
            if let Some(chirality) = chirality {
                builder = builder.with_chirality(chirality);
            }
            builder.build()
        }
    }
}

fn atom_with_hydrogen_count(atom: Atom, hydrogens: u8) -> Atom {
    match atom.syntax() {
        AtomSyntax::OrganicSubset if hydrogens == atom.hydrogen_count() => atom,
        AtomSyntax::OrganicSubset | AtomSyntax::Bracket => {
            let mut builder = Atom::builder()
                .with_symbol(atom.symbol())
                .with_aromatic(atom.aromatic())
                .with_hydrogens(hydrogens)
                .with_charge(atom.charge())
                .with_class(atom.class());
            if let Some(isotope) = atom.isotope_mass_number() {
                builder = builder.with_isotope(isotope);
            }
            if let Some(chirality) = atom.chirality() {
                builder = builder.with_chirality(chirality);
            }
            builder.build()
        }
    }
}

fn maybe_collapse_atom_to_organic_subset(smiles: &Smiles, node_id: usize, atom: Atom) -> Atom {
    if atom.syntax() != AtomSyntax::Bracket
        || atom.isotope_mass_number().is_some()
        || atom.charge_value() != 0
        || atom.class() != 0
        || atom.chirality().is_some()
        || !canonicalization_valid_unbracketed(atom.symbol())
        || implicit_hydrogens_if_written_unbracketed(smiles, node_id, &atom)
            != atom.hydrogen_count()
    {
        return atom;
    }

    Atom::new_organic_subset(atom.symbol(), atom.aromatic())
}

fn canonicalization_valid_unbracketed(symbol: AtomSymbol) -> bool {
    matches!(
        symbol,
        AtomSymbol::WildCard
            | AtomSymbol::Element(
                Element::B
                    | Element::C
                    | Element::N
                    | Element::O
                    | Element::P
                    | Element::S
                    | Element::F
                    | Element::Cl
                    | Element::Br
                    | Element::I
            )
    )
}

fn canonical_stereo_neighbor_sort_key(
    neighbor: StereoNeighbor,
    new_index_of_old_node: &[usize],
) -> (u8, usize) {
    match neighbor {
        StereoNeighbor::Atom(node_id) => (0, new_index_of_old_node[node_id]),
        StereoNeighbor::ExplicitHydrogen => (1, usize::MAX),
    }
}

fn allene_like_stereo_center(smiles: &Smiles, node_id: usize) -> bool {
    let incident_bonds = smiles.edges_for_node(node_id);
    incident_bonds.len() == 2 && incident_bonds.iter().all(|edge| edge.bond() == Bond::Double)
}

fn non_tetrahedral_default_chirality(
    chirality: Option<Chirality>,
    parsed_neighbors: &[StereoNeighbor],
) -> Option<Chirality> {
    match chirality {
        Some(Chirality::At) if parsed_neighbors.len() == 5 => Some(Chirality::TB(1)),
        Some(Chirality::AtAt) if parsed_neighbors.len() == 5 => Some(Chirality::TB(2)),
        Some(Chirality::At) if parsed_neighbors.len() == 6 => Some(Chirality::OH(1)),
        Some(Chirality::AtAt) if parsed_neighbors.len() == 6 => Some(Chirality::OH(2)),
        other => other,
    }
}

fn stereo_chirality_normal_form(
    smiles: &Smiles,
    node_id: usize,
    chirality: Option<Chirality>,
    parsed_neighbors: &[StereoNeighbor],
) -> Option<Chirality> {
    match chirality {
        Some(Chirality::At) => {
            if allene_like_stereo_center(smiles, node_id) {
                Some(Chirality::AL(1))
            } else if parsed_neighbors.len() == 4 {
                Some(Chirality::TH(1))
            } else if parsed_neighbors.len() == 5 {
                Some(Chirality::TB(1))
            } else if parsed_neighbors.len() == 6 {
                Some(Chirality::OH(1))
            } else {
                Some(Chirality::At)
            }
        }
        Some(Chirality::AtAt) => {
            if allene_like_stereo_center(smiles, node_id) {
                Some(Chirality::AL(2))
            } else if parsed_neighbors.len() == 4 {
                Some(Chirality::TH(2))
            } else if parsed_neighbors.len() == 5 {
                Some(Chirality::TB(2))
            } else if parsed_neighbors.len() == 6 {
                Some(Chirality::OH(2))
            } else {
                Some(Chirality::AtAt)
            }
        }
        other => non_tetrahedral_default_chirality(other, parsed_neighbors),
    }
}

fn supports_stereo_neighbor_normalization(chirality: Option<Chirality>) -> bool {
    matches!(
        chirality,
        Some(
            Chirality::At
                | Chirality::AtAt
                | Chirality::TH(_)
                | Chirality::AL(_)
                | Chirality::SP(_)
                | Chirality::TB(_)
                | Chirality::OH(_)
        )
    )
}

fn expected_stereo_neighbor_count(chirality: Chirality) -> usize {
    match chirality {
        Chirality::At
        | Chirality::AtAt
        | Chirality::TH(_)
        | Chirality::AL(_)
        | Chirality::SP(_) => 4,
        Chirality::TB(_) => 5,
        Chirality::OH(_) => 6,
    }
}

fn stereo_neighbors_with_implicit_hydrogens(
    smiles: &Smiles,
    node_id: usize,
    chirality: Chirality,
    parsed_neighbors: &[StereoNeighbor],
) -> Vec<StereoNeighbor> {
    let mut neighbors = parsed_neighbors.to_vec();
    let expected_count = expected_stereo_neighbor_count(chirality);
    let atom_hydrogens = usize::from(smiles.nodes()[node_id].hydrogen_count());
    let present_implicit_hydrogens =
        neighbors.iter().filter(|neighbor| **neighbor == StereoNeighbor::ExplicitHydrogen).count();
    let missing_implicit_hydrogens = atom_hydrogens.saturating_sub(present_implicit_hydrogens);
    let missing_capacity = expected_count.saturating_sub(neighbors.len());

    for _ in 0..missing_implicit_hydrogens.min(missing_capacity) {
        neighbors.push(StereoNeighbor::ExplicitHydrogen);
    }
    neighbors
}

fn canonical_stereo_neighbors_row(
    smiles: &Smiles,
    node_id: usize,
    chirality: Option<Chirality>,
    parsed_neighbors: &[StereoNeighbor],
    new_index_of_old_node: &[usize],
) -> Vec<StereoNeighbor> {
    let chirality = stereo_chirality_normal_form(smiles, node_id, chirality, parsed_neighbors);
    if !supports_stereo_neighbor_normalization(chirality) {
        return parsed_neighbors.to_vec();
    }
    let expanded_neighbors = stereo_neighbors_with_implicit_hydrogens(
        smiles,
        node_id,
        chirality.unwrap_or_else(|| unreachable!()),
        parsed_neighbors,
    );

    match chirality {
        Some(Chirality::At | Chirality::AtAt | Chirality::TH(_) | Chirality::AL(_)) => {
            let mut neighbors = expanded_neighbors;
            neighbors.sort_unstable_by_key(|neighbor| {
                canonical_stereo_neighbor_sort_key(*neighbor, new_index_of_old_node)
            });
            neighbors
        }
        Some(Chirality::SP(kind)) => {
            square_assignment_from_shape(Chirality::SP(kind), &expanded_neighbors).map_or_else(
                || parsed_neighbors.to_vec(),
                |assignment| canonical_square_planar_neighbors(&assignment, new_index_of_old_node),
            )
        }
        Some(Chirality::TB(_)) => {
            canonical_tb_neighbors(chirality, &expanded_neighbors, new_index_of_old_node)
                .unwrap_or_else(|| parsed_neighbors.to_vec())
        }
        Some(Chirality::OH(_)) => {
            canonical_oh_neighbors(chirality, &expanded_neighbors, new_index_of_old_node)
                .unwrap_or_else(|| parsed_neighbors.to_vec())
        }
        None => parsed_neighbors.to_vec(),
    }
}

fn normalized_stereo_chirality(
    smiles: &Smiles,
    node_id: usize,
    chirality: Option<Chirality>,
    parsed_neighbors: &[StereoNeighbor],
    normalized_neighbors: &[StereoNeighbor],
    rooted_classes: &[usize],
    refined_classes: &[usize],
) -> Option<Chirality> {
    let chirality = stereo_chirality_normal_form(smiles, node_id, chirality, parsed_neighbors)?;
    match chirality {
        Chirality::At | Chirality::AtAt | Chirality::TH(_) | Chirality::AL(_) => {
            if matches!(chirality, Chirality::At | Chirality::AtAt | Chirality::TH(_))
                && !tetrahedral_like_is_stereogenic(
                    smiles,
                    node_id,
                    chirality,
                    parsed_neighbors,
                    rooted_classes,
                    refined_classes,
                )
            {
                return None;
            }
            normalized_tetrahedral_chirality(
                Some(chirality),
                parsed_neighbors,
                normalized_neighbors,
            )
        }
        Chirality::SP(_) => Some(Chirality::SP(1)),
        Chirality::TB(_) => Some(Chirality::TB(1)),
        Chirality::OH(_) => Some(Chirality::OH(1)),
    }
}

fn stereo_substituent_identity_key(
    smiles: &Smiles,
    endpoint: usize,
    neighbor: StereoNeighbor,
    rooted_classes: &[usize],
    refined_classes: &[usize],
) -> StereoSubstituentIdentityKey {
    match neighbor {
        StereoNeighbor::ExplicitHydrogen => StereoSubstituentIdentityKey::ExplicitHydrogen,
        StereoNeighbor::Atom(node_id) => {
            StereoSubstituentIdentityKey::Atom(atom_based_substituent_priority_key(
                smiles,
                endpoint,
                node_id,
                rooted_classes,
                refined_classes,
            ))
        }
    }
}

fn tetrahedral_like_is_stereogenic(
    smiles: &Smiles,
    node_id: usize,
    chirality: Chirality,
    parsed_neighbors: &[StereoNeighbor],
    rooted_classes: &[usize],
    refined_classes: &[usize],
) -> bool {
    let mut keys =
        stereo_neighbors_with_implicit_hydrogens(smiles, node_id, chirality, parsed_neighbors)
            .into_iter()
            .map(|neighbor| {
                stereo_substituent_identity_key(
                    smiles,
                    node_id,
                    neighbor,
                    rooted_classes,
                    refined_classes,
                )
            })
            .collect::<Vec<_>>();
    keys.sort_unstable();
    keys.windows(2).all(|window| window[0] != window[1])
}

fn canonicalization_atom_spelling_normal_form(smiles: &Smiles, node_id: usize, atom: Atom) -> Atom {
    if atom.syntax() == AtomSyntax::Bracket {
        maybe_collapse_atom_to_organic_subset(smiles, node_id, atom)
    } else {
        atom
    }
}

fn canonical_stereo_neighbor_key(neighbor: StereoNeighbor) -> CanonicalStereoNeighborKey {
    match neighbor {
        StereoNeighbor::Atom(node_id) => CanonicalStereoNeighborKey(0, node_id),
        StereoNeighbor::ExplicitHydrogen => CanonicalStereoNeighborKey(1, usize::MAX),
    }
}

fn stereo_neighbor_sequence_key(
    neighbors: &[StereoNeighbor],
    new_index_of_old_node: &[usize],
) -> Vec<(u8, usize)> {
    neighbors
        .iter()
        .copied()
        .map(|neighbor| canonical_stereo_neighbor_sort_key(neighbor, new_index_of_old_node))
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SquareShape {
    U,
    Z,
    Four,
}

const fn square_shape_path(shape: SquareShape) -> [usize; 4] {
    match shape {
        SquareShape::U => [0, 1, 2, 3],
        SquareShape::Z => [0, 1, 3, 2],
        SquareShape::Four => [0, 2, 1, 3],
    }
}

fn square_shape_for_chirality(chirality: Chirality) -> Option<SquareShape> {
    match chirality {
        Chirality::SP(1) => Some(SquareShape::U),
        Chirality::SP(2) => Some(SquareShape::Four),
        Chirality::SP(3) => Some(SquareShape::Z),
        _ => None,
    }
}

fn canonical_square_planar_neighbors(
    parsed_neighbors: &[StereoNeighbor],
    new_index_of_old_node: &[usize],
) -> Vec<StereoNeighbor> {
    let assignment =
        [parsed_neighbors[0], parsed_neighbors[1], parsed_neighbors[2], parsed_neighbors[3]];
    let mut best = assignment;
    let mut best_key = stereo_neighbor_sequence_key(&best, new_index_of_old_node);
    for rotation in 0..4 {
        let rotated = [
            assignment[(rotation) % 4],
            assignment[(rotation + 1) % 4],
            assignment[(rotation + 2) % 4],
            assignment[(rotation + 3) % 4],
        ];
        let rotated_key = stereo_neighbor_sequence_key(&rotated, new_index_of_old_node);
        if rotated_key < best_key {
            best = rotated;
            best_key = rotated_key;
        }

        let mirrored = [
            assignment[(rotation) % 4],
            assignment[(rotation + 3) % 4],
            assignment[(rotation + 2) % 4],
            assignment[(rotation + 1) % 4],
        ];
        let mirrored_key = stereo_neighbor_sequence_key(&mirrored, new_index_of_old_node);
        if mirrored_key < best_key {
            best = mirrored;
            best_key = mirrored_key;
        }
    }
    best.to_vec()
}

fn square_assignment_from_shape(
    chirality: Chirality,
    sequence: &[StereoNeighbor],
) -> Option<[StereoNeighbor; 4]> {
    let path = square_shape_path(square_shape_for_chirality(chirality)?);
    let mut assignment = [StereoNeighbor::ExplicitHydrogen; 4];
    for (index, &position) in path.iter().enumerate() {
        assignment[position] =
            sequence.get(index).copied().unwrap_or(StereoNeighbor::ExplicitHydrogen);
    }
    Some(assignment)
}

fn tb_axis_and_order(chirality: Chirality) -> Option<(usize, usize, bool)> {
    match chirality {
        Chirality::TB(1) => Some((0, 4, false)),
        Chirality::TB(2) => Some((0, 4, true)),
        Chirality::TB(3) => Some((0, 3, false)),
        Chirality::TB(4) => Some((0, 3, true)),
        Chirality::TB(5) => Some((0, 2, false)),
        Chirality::TB(6) => Some((0, 2, true)),
        Chirality::TB(7) => Some((0, 1, false)),
        Chirality::TB(8) => Some((0, 1, true)),
        Chirality::TB(9) => Some((1, 4, false)),
        Chirality::TB(11) => Some((1, 4, true)),
        Chirality::TB(10) => Some((1, 3, false)),
        Chirality::TB(12) => Some((1, 3, true)),
        Chirality::TB(13) => Some((1, 2, false)),
        Chirality::TB(14) => Some((1, 2, true)),
        Chirality::TB(15) => Some((2, 4, false)),
        Chirality::TB(20) => Some((2, 4, true)),
        Chirality::TB(16) => Some((2, 3, false)),
        Chirality::TB(19) => Some((2, 3, true)),
        Chirality::TB(17) => Some((3, 4, false)),
        Chirality::TB(18) => Some((3, 4, true)),
        _ => None,
    }
}

fn canonical_tb_neighbors(
    chirality: Option<Chirality>,
    parsed_neighbors: &[StereoNeighbor],
    new_index_of_old_node: &[usize],
) -> Option<Vec<StereoNeighbor>> {
    let chirality = chirality?;
    let (axis_start_index, axis_end_index, clockwise) = tb_axis_and_order(chirality)?;
    let axis_start = parsed_neighbors.get(axis_start_index).copied()?;
    let axis_end = parsed_neighbors.get(axis_end_index).copied()?;
    let mut cycle = parsed_neighbors
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(index, neighbor)| {
            (index != axis_start_index && index != axis_end_index).then_some(neighbor)
        })
        .collect::<Vec<_>>();
    if clockwise {
        cycle.reverse();
    }
    if canonical_stereo_neighbor_sort_key(axis_end, new_index_of_old_node)
        < canonical_stereo_neighbor_sort_key(axis_start, new_index_of_old_node)
    {
        cycle.reverse();
        return Some(
            [
                vec![axis_end],
                minimal_cycle_rotation(&cycle, new_index_of_old_node),
                vec![axis_start],
            ]
            .concat(),
        );
    }
    Some(
        [vec![axis_start], minimal_cycle_rotation(&cycle, new_index_of_old_node), vec![axis_end]]
            .concat(),
    )
}

fn canonical_oh_neighbors(
    chirality: Option<Chirality>,
    parsed_neighbors: &[StereoNeighbor],
    new_index_of_old_node: &[usize],
) -> Option<Vec<StereoNeighbor>> {
    let permutation = octahedral_normalization_permutation(chirality?)?;
    let normalized = permutation
        .iter()
        .copied()
        .map(|index| parsed_neighbors.get(index).copied())
        .collect::<Option<Vec<_>>>()?;
    let mut best = normalized.clone();
    let mut best_key = stereo_neighbor_sequence_key(&best, new_index_of_old_node);

    for permutation in OCTAHEDRAL_OH1_EQUIVALENT_PERMUTATIONS {
        let candidate =
            permutation.iter().copied().map(|index| normalized[index]).collect::<Vec<_>>();
        let candidate_key = stereo_neighbor_sequence_key(&candidate, new_index_of_old_node);
        if candidate_key < best_key {
            best = candidate;
            best_key = candidate_key;
        }
    }

    Some(best)
}

fn minimal_cycle_rotation(
    cycle: &[StereoNeighbor],
    new_index_of_old_node: &[usize],
) -> Vec<StereoNeighbor> {
    let mut best = cycle.to_vec();
    let mut best_key = stereo_neighbor_sequence_key(&best, new_index_of_old_node);
    for rotation in 1..cycle.len() {
        let rotated =
            cycle[rotation..].iter().chain(cycle[..rotation].iter()).copied().collect::<Vec<_>>();
        let rotated_key = stereo_neighbor_sequence_key(&rotated, new_index_of_old_node);
        if rotated_key < best_key {
            best = rotated;
            best_key = rotated_key;
        }
    }
    best
}

const fn octahedral_normalization_permutation(chirality: Chirality) -> Option<[usize; 6]> {
    match chirality {
        Chirality::OH(1) => Some([0, 1, 2, 3, 4, 5]),
        Chirality::OH(2) => Some([0, 1, 4, 3, 2, 5]),
        Chirality::OH(3) => Some([0, 1, 2, 3, 5, 4]),
        Chirality::OH(4) => Some([0, 1, 2, 4, 3, 5]),
        Chirality::OH(5) => Some([0, 1, 2, 5, 3, 4]),
        Chirality::OH(6) => Some([0, 1, 2, 4, 5, 3]),
        Chirality::OH(7) => Some([0, 1, 2, 5, 4, 3]),
        Chirality::OH(8) => Some([0, 1, 3, 2, 4, 5]),
        Chirality::OH(9) => Some([0, 1, 3, 2, 5, 4]),
        Chirality::OH(10) => Some([0, 1, 4, 2, 3, 5]),
        Chirality::OH(11) => Some([0, 1, 5, 2, 3, 4]),
        Chirality::OH(12) => Some([0, 1, 4, 2, 5, 3]),
        Chirality::OH(13) => Some([0, 1, 5, 2, 4, 3]),
        Chirality::OH(14) => Some([0, 1, 3, 4, 2, 5]),
        Chirality::OH(15) => Some([0, 1, 3, 5, 2, 4]),
        Chirality::OH(16) => Some([0, 1, 5, 3, 2, 4]),
        Chirality::OH(17) => Some([0, 1, 4, 5, 2, 3]),
        Chirality::OH(18) => Some([0, 1, 5, 4, 2, 3]),
        Chirality::OH(19) => Some([0, 1, 3, 4, 5, 2]),
        Chirality::OH(20) => Some([0, 1, 3, 5, 4, 2]),
        Chirality::OH(21) => Some([0, 1, 4, 3, 5, 2]),
        Chirality::OH(22) => Some([0, 1, 5, 3, 4, 2]),
        Chirality::OH(23) => Some([0, 1, 4, 5, 3, 2]),
        Chirality::OH(24) => Some([0, 1, 5, 4, 3, 2]),
        Chirality::OH(25) => Some([0, 2, 3, 4, 5, 1]),
        Chirality::OH(26) => Some([0, 2, 3, 5, 4, 1]),
        Chirality::OH(27) => Some([0, 2, 4, 3, 5, 1]),
        Chirality::OH(28) => Some([0, 2, 5, 3, 4, 1]),
        Chirality::OH(29) => Some([0, 2, 4, 5, 3, 1]),
        Chirality::OH(30) => Some([0, 2, 5, 4, 3, 1]),
        _ => None,
    }
}

const OCTAHEDRAL_OH1_EQUIVALENT_PERMUTATIONS: [[usize; 6]; 24] = [
    [0, 1, 2, 3, 4, 5],
    [0, 2, 3, 4, 1, 5],
    [0, 3, 4, 1, 2, 5],
    [0, 4, 1, 2, 3, 5],
    [1, 0, 4, 5, 2, 3],
    [1, 2, 0, 4, 5, 3],
    [1, 4, 5, 2, 0, 3],
    [1, 5, 2, 0, 4, 3],
    [2, 0, 1, 5, 3, 4],
    [2, 1, 5, 3, 0, 4],
    [2, 3, 0, 1, 5, 4],
    [2, 5, 3, 0, 1, 4],
    [3, 0, 2, 5, 4, 1],
    [3, 2, 5, 4, 0, 1],
    [3, 4, 0, 2, 5, 1],
    [3, 5, 4, 0, 2, 1],
    [4, 0, 3, 5, 1, 2],
    [4, 1, 0, 3, 5, 2],
    [4, 3, 5, 1, 0, 2],
    [4, 5, 1, 0, 3, 2],
    [5, 1, 4, 3, 2, 0],
    [5, 2, 1, 4, 3, 0],
    [5, 3, 2, 1, 4, 0],
    [5, 4, 3, 2, 1, 0],
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct AtomBasedSubstituentPriorityKey {
    atomic_number: u8,
    isotope_mass_number: Option<u16>,
    charge: i8,
    aromatic: bool,
    bond_order_to_endpoint: u8,
    rooted_class: usize,
    refined_class: usize,
}

fn atom_based_substituent_priority_key(
    smiles: &Smiles,
    endpoint: usize,
    neighbor: usize,
    rooted_classes: &[usize],
    refined_classes: &[usize],
) -> AtomBasedSubstituentPriorityKey {
    let atom = smiles.node_by_id(neighbor).unwrap_or_else(|| unreachable!());
    let atomic_number = atom.element().map_or(0, u8::from);
    let bond_order_to_endpoint = match smiles
        .edge_for_node_pair((endpoint, neighbor))
        .unwrap_or_else(|| unreachable!())
        .bond()
    {
        Bond::Single | Bond::Up | Bond::Down | Bond::Aromatic => 1,
        Bond::Double => 2,
        Bond::Triple => 3,
        Bond::Quadruple => 4,
    };

    AtomBasedSubstituentPriorityKey {
        atomic_number,
        isotope_mass_number: atom.isotope_mass_number(),
        charge: atom.charge_value(),
        aromatic: atom.aromatic(),
        bond_order_to_endpoint,
        rooted_class: rooted_classes[neighbor],
        refined_class: refined_classes[neighbor],
    }
}

fn non_semantic_double_bond_supports_semantic_stereo(
    smiles: &Smiles,
    node_a: usize,
    node_b: usize,
) -> bool {
    smiles
        .bond_matrix()
        .sparse_row_values_ref(node_a)
        .filter(|entry| !matches!(entry.bond(), Bond::Single | Bond::Up | Bond::Down))
        .count()
        == 1
        && smiles
            .bond_matrix()
            .sparse_row_values_ref(node_b)
            .filter(|entry| !matches!(entry.bond(), Bond::Single | Bond::Up | Bond::Down))
            .count()
            == 1
}

fn non_semantic_double_bond_is_in_cycle(smiles: &Smiles, node_a: usize, node_b: usize) -> bool {
    let mut queue = alloc::collections::VecDeque::from([node_a]);
    let mut seen = vec![false; smiles.nodes().len()];
    seen[node_a] = true;

    while let Some(current) = queue.pop_front() {
        for neighbor in smiles.bond_matrix().sparse_row(current) {
            if (current == node_a && neighbor == node_b)
                || (current == node_b && neighbor == node_a)
            {
                continue;
            }
            if neighbor == node_b {
                return true;
            }
            if !seen[neighbor] {
                seen[neighbor] = true;
                queue.push_back(neighbor);
            }
        }
    }

    false
}

fn non_semantic_double_bond_directional_neighbors(
    smiles: &Smiles,
    endpoint: usize,
    opposite_endpoint: usize,
) -> Vec<usize> {
    smiles
        .bond_matrix()
        .sparse_row(endpoint)
        .zip(smiles.bond_matrix().sparse_row_values_ref(endpoint))
        .filter_map(|(neighbor, entry)| {
            (neighbor != opposite_endpoint && matches!(entry.bond(), Bond::Up | Bond::Down))
                .then_some(neighbor)
        })
        .collect()
}

fn non_semantic_double_bond_has_unique_reference_substituent(
    smiles: &Smiles,
    endpoint: usize,
    opposite_endpoint: usize,
    rooted_classes: &[usize],
    refined_classes: &[usize],
) -> Option<usize> {
    let neighbors = smiles
        .bond_matrix()
        .sparse_row(endpoint)
        .filter(|&neighbor| neighbor != opposite_endpoint)
        .collect::<Vec<_>>();
    let (&first, rest) = neighbors.split_first()?;
    let mut best = first;
    let mut best_key = atom_based_substituent_priority_key(
        smiles,
        endpoint,
        best,
        rooted_classes,
        refined_classes,
    );
    let mut unique_best = true;

    for &candidate in rest {
        let candidate_key = atom_based_substituent_priority_key(
            smiles,
            endpoint,
            candidate,
            rooted_classes,
            refined_classes,
        );
        match candidate_key.cmp(&best_key) {
            core::cmp::Ordering::Greater => {
                best = candidate;
                best_key = candidate_key;
                unique_best = true;
            }
            core::cmp::Ordering::Equal => unique_best = false,
            core::cmp::Ordering::Less => {}
        }
    }

    unique_best.then_some(best)
}

fn atom_based_chirality_is_clockwise(chirality: Chirality) -> bool {
    matches!(chirality, Chirality::AtAt | Chirality::TH(2))
}

fn atom_based_double_bond_target_neighbors(
    parsed_neighbors: &[StereoNeighbor],
    opposite_endpoint: usize,
    reference_atom: usize,
) -> Option<Vec<StereoNeighbor>> {
    let reference = StereoNeighbor::Atom(reference_atom);
    let opposite = StereoNeighbor::Atom(opposite_endpoint);
    let other = parsed_neighbors
        .iter()
        .copied()
        .find(|neighbor| *neighbor != reference && *neighbor != opposite)?;
    Some(vec![reference, opposite, other])
}

fn atom_based_permutation_from(
    parsed_neighbors: &[StereoNeighbor],
    target_neighbors: &[StereoNeighbor],
) -> Option<Vec<usize>> {
    let mut used = vec![false; parsed_neighbors.len()];
    let mut permutation = Vec::with_capacity(target_neighbors.len());
    for target in target_neighbors {
        let index = parsed_neighbors
            .iter()
            .enumerate()
            .find_map(|(index, parsed)| (!used[index] && parsed == target).then_some(index))?;
        used[index] = true;
        permutation.push(index);
    }
    Some(permutation)
}

fn atom_based_permutation_is_odd(permutation: &[usize]) -> bool {
    let mut visited = vec![false; permutation.len()];
    let mut transpositions = 0;

    for start in 0..permutation.len() {
        if visited[start] {
            continue;
        }
        let mut length = 0;
        let mut cursor = start;
        while !visited[cursor] {
            visited[cursor] = true;
            cursor = permutation[cursor];
            length += 1;
        }
        if length > 0 {
            transpositions += length - 1;
        }
    }

    transpositions % 2 == 1
}

fn canonical_directional_edge_key(
    edge_key: (usize, usize),
    preorder_indices: &[usize],
) -> (usize, usize) {
    let left = preorder_indices[edge_key.0];
    let right = preorder_indices[edge_key.1];
    if left <= right { (left, right) } else { (right, left) }
}

fn atom_based_override_bond(rows: &[Vec<(usize, Bond)>], from: usize, to: usize) -> Option<Bond> {
    let (left, right) = Smiles::edge_key(from, to);
    let row = rows.get(left)?;
    let index = row.binary_search_by_key(&right, |&(neighbor, _)| neighbor).ok()?;
    Some(row[index].1)
}

fn canonicalization_state_key(smiles: &Smiles) -> CanonicalizationStateKey {
    let atom_labels = smiles.nodes().iter().copied().map(canonical_atom_label).collect();
    let bond_edges = smiles
        .bond_matrix()
        .sparse_entries()
        .filter(|((row, column), _)| row < column)
        .map(|((row, column), entry)| (row, column, canonical_bond_label(entry.bond())))
        .collect();
    let parsed_stereo_neighbors = smiles
        .parsed_stereo_neighbors
        .iter()
        .map(|row| row.iter().copied().map(canonical_stereo_neighbor_key).collect())
        .collect();

    CanonicalizationStateKey {
        atom_labels,
        bond_edges,
        parsed_stereo_neighbors,
        implicit_hydrogen_cache: smiles.implicit_hydrogen_cache.clone(),
    }
}

const fn canonical_bond_label(bond: Bond) -> CanonicalBondLabel {
    CanonicalBondLabel(match bond {
        Bond::Single => 0,
        Bond::Double => 1,
        Bond::Triple => 2,
        Bond::Quadruple => 3,
        Bond::Aromatic => 4,
        Bond::Up => 5,
        Bond::Down => 6,
    })
}

const fn canonical_chirality_key(chirality: Option<Chirality>) -> (u8, u8) {
    match chirality {
        None => (0, 0),
        Some(Chirality::At) => (1, 0),
        Some(Chirality::AtAt) => (2, 0),
        Some(Chirality::TH(value)) => (3, value),
        Some(Chirality::AL(value)) => (4, value),
        Some(Chirality::SP(value)) => (5, value),
        Some(Chirality::TB(value)) => (6, value),
        Some(Chirality::OH(value)) => (7, value),
    }
}

fn remap_parsed_stereo_neighbors_row(
    smiles: &Smiles,
    old_node: usize,
    new_index_of_old_node: &[usize],
) -> Vec<StereoNeighbor> {
    smiles
        .parsed_stereo_neighbors
        .get(old_node)
        .into_iter()
        .flatten()
        .copied()
        .map(|neighbor| {
            match neighbor {
                StereoNeighbor::Atom(old_neighbor) => {
                    StereoNeighbor::Atom(new_index_of_old_node[old_neighbor])
                }
                StereoNeighbor::ExplicitHydrogen => StereoNeighbor::ExplicitHydrogen,
            }
        })
        .collect()
}
