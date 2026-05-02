//! Module for bonds stored as graph edge tuples.

use crate::bond::{Bond, BondDescriptor, ring_num::RingNum};

/// Contains the two atom indices connected via the [`Bond`].
pub type BondEdge = (usize, usize, Bond, Option<RingNum>, bool);

/// Creates a new edge tuple.
///
/// # Examples
///
/// ```
/// use smiles_parser::bond::{Bond, bond_edge::bond_edge};
///
/// let edge = bond_edge(0, 1, Bond::Double, None);
/// assert_eq!(edge, (0, 1, Bond::Double, None, false));
/// ```
#[inline]
#[must_use]
pub const fn bond_edge(
    node_a: usize,
    node_b: usize,
    bond: Bond,
    ring_num: Option<RingNum>,
) -> BondEdge {
    bond_edge_with_aromaticity(node_a, node_b, bond, ring_num, false)
}

/// Creates a new edge tuple with explicit aromaticity.
#[inline]
#[must_use]
pub const fn bond_edge_with_aromaticity(
    node_a: usize,
    node_b: usize,
    bond: Bond,
    ring_num: Option<RingNum>,
    aromatic: bool,
) -> BondEdge {
    (node_a, node_b, bond, ring_num, aromatic)
}

/// Creates a new edge tuple from a bond descriptor.
#[inline]
#[must_use]
pub const fn bond_edge_from_descriptor(
    node_a: usize,
    node_b: usize,
    descriptor: BondDescriptor,
    ring_num: Option<RingNum>,
) -> BondEdge {
    bond_edge_with_aromaticity(
        node_a,
        node_b,
        descriptor.bond(),
        ring_num,
        descriptor.is_aromatic(),
    )
}

/// Returns the other node id for the provided incident edge, if any.
///
/// # Examples
///
/// ```
/// use smiles_parser::bond::{
///     Bond,
///     bond_edge::{bond_edge, bond_edge_other},
/// };
///
/// let edge = bond_edge(2, 5, Bond::Single, None);
/// assert_eq!(bond_edge_other(edge, 2), Some(5));
/// assert_eq!(bond_edge_other(edge, 99), None);
/// ```
#[inline]
#[must_use]
pub const fn bond_edge_other(edge: BondEdge, node_id: usize) -> Option<usize> {
    if edge.0 == node_id {
        Some(edge.1)
    } else if edge.1 == node_id {
        Some(edge.0)
    } else {
        None
    }
}

/// Returns the ring number value stored in the edge, if any.
///
/// # Examples
///
/// ```
/// use smiles_parser::bond::{
///     Bond,
///     bond_edge::{bond_edge, bond_edge_ring_num_val},
///     ring_num::RingNum,
/// };
///
/// let edge = bond_edge(0, 1, Bond::Single, Some(RingNum::try_new(7)?));
/// assert_eq!(bond_edge_ring_num_val(edge), Some(7));
/// # Ok::<(), smiles_parser::SmilesError>(())
/// ```
#[inline]
#[must_use]
pub fn bond_edge_ring_num_val(edge: BondEdge) -> Option<u8> {
    edge.3.map(|num| num.get())
}

#[cfg(test)]
mod tests {
    use crate::bond::{
        Bond, BondDescriptor,
        bond_edge::{
            BondEdge, bond_edge, bond_edge_from_descriptor, bond_edge_other,
            bond_edge_ring_num_val, bond_edge_with_aromaticity,
        },
        ring_num::RingNum,
    };

    #[test]
    fn test_bond_edge_new_and_accessors() {
        let edge = bond_edge(3, 7, Bond::Double, None);

        assert_eq!(edge.0, 3);
        assert_eq!(edge.1, 7);
        assert_eq!((edge.0, edge.1), (3, 7));
        assert_eq!(edge.2, Bond::Double);
        assert_eq!(edge.3, None);
        assert!(!edge.4);

        let aromatic = bond_edge_with_aromaticity(3, 7, Bond::Triple, None, true);
        assert_eq!(aromatic.2, Bond::Triple);
        assert!(aromatic.4);

        let from_descriptor =
            bond_edge_from_descriptor(3, 7, BondDescriptor::aromatic(Bond::Single), None);
        assert_eq!(from_descriptor, bond_edge_with_aromaticity(3, 7, Bond::Single, None, true));
    }

    #[test]
    fn test_bond_edge_other() {
        let edge = bond_edge(10, 2, Bond::Single, None);
        assert_eq!(bond_edge_other(edge, 10), Some(2));
        assert_eq!(bond_edge_other(edge, 2), Some(10));
        assert_eq!(bond_edge_other(edge, 5), None);
    }

    #[test]
    fn test_bond_edge_ring_num_value() {
        let ring_num = RingNum::try_new(2).unwrap();
        let edge = bond_edge(0, 1, Bond::Single, Some(ring_num));
        assert_eq!(bond_edge_ring_num_val(edge), Some(2));
        assert_eq!(bond_edge_ring_num_val((0, 1, Bond::Single, None, false)), None);
    }

    #[test]
    fn test_bond_edge_vertices_preserve_order() {
        let edge: BondEdge = bond_edge(10, 2, Bond::Single, None);
        assert_eq!((edge.0, edge.1), (10, 2));
    }
}
