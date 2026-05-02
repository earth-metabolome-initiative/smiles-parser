#![no_std]
#![doc = include_str!("../README.md")]

#[allow(unused_imports)]
#[macro_use]
extern crate alloc;
#[cfg(test)]
#[macro_use]
extern crate std;
#[cfg(all(feature = "datasets", not(test)))]
extern crate std;

pub mod atom;
pub mod bond;
#[cfg(feature = "datasets")]
pub mod datasets;
pub mod errors;
pub(crate) mod parser;
pub mod smiles;
pub mod token;

#[cfg(feature = "datasets")]
pub use crate::datasets::{
    CacheMode, DatasetArtifact, DatasetCompression, DatasetError, DatasetFetchOptions,
    DatasetSmilesIter, DatasetSource, GzipMode, MASS_SPEC_GYM_SMILES, MassSpecGymSmiles,
    PUBCHEM_SMILES, PubChemSmiles, SmilesDatasetSource, default_dataset_cache_dir,
};
pub use crate::{
    errors::{SmilesError, SmilesErrorWithSpan},
    smiles::{
        AromaticityAssignment, AromaticityAssignmentApplicationError, AromaticityDiagnostic,
        AromaticityModel, AromaticityPerception, AromaticityPolicy, AromaticityRingFamilyKind,
        AromaticityStatus, DoubleBondStereoConfig, KekulizationError, KekulizationMode,
        RdkitDefaultAromaticity, RdkitMdlAromaticity, RdkitSimpleAromaticity, RingAtomMembership,
        RingAtomMembershipScratch, RingMembership, Smiles, SmilesComponents, SymmSssrResult,
        SymmSssrStatus, WildcardAromaticityPerception, WildcardMolecularFormulaConversionError,
        WildcardSmiles, WildcardSmilesComponents,
    },
};

/// Common imports for working with this crate.
pub mod prelude {
    pub use crate::{
        AromaticityAssignment, AromaticityAssignmentApplicationError, AromaticityDiagnostic,
        AromaticityModel, AromaticityPerception, AromaticityPolicy, AromaticityRingFamilyKind,
        AromaticityStatus, DoubleBondStereoConfig, KekulizationError, KekulizationMode,
        RdkitDefaultAromaticity, RdkitMdlAromaticity, RdkitSimpleAromaticity, RingAtomMembership,
        RingAtomMembershipScratch, RingMembership, Smiles, SmilesComponents, SmilesError,
        SmilesErrorWithSpan, SymmSssrResult, SymmSssrStatus, WildcardAromaticityPerception,
        WildcardMolecularFormulaConversionError, WildcardSmiles, WildcardSmilesComponents,
    };
    #[cfg(feature = "datasets")]
    pub use crate::{
        CacheMode, DatasetArtifact, DatasetCompression, DatasetError, DatasetFetchOptions,
        DatasetSmilesIter, DatasetSource, GzipMode, MASS_SPEC_GYM_SMILES, MassSpecGymSmiles,
        PUBCHEM_SMILES, PubChemSmiles, SmilesDatasetSource, default_dataset_cache_dir,
    };
}
