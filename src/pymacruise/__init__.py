"""Python-native MaCRUISE workflow."""

from pymacruise.image_io import ImageGeometryError, load_nifti_ras, save_nifti_ras
from pymacruise.raw_preparation import RawPreparationResult, prepare_raw
from pymacruise.contracts import (
    ValidatedInputs,
    WorkflowInputError,
    WorkflowProfile,
    validate_workflow_files,
    validate_workflow_inputs,
)
from pymacruise.cortical_separation import cortical_separation
from pymacruise.boundary import (
    clean_csf_boundary,
    clean_csf_boundary_image,
    combine_multi_atlas_toads,
    combine_multi_atlas_toads_images,
    merge_csf_boundary,
    merge_csf_boundary_images,
    multi_atlas_csf_boundary,
    multi_atlas_csf_boundary_image,
)
from pymacruise.cropping import CropParameters, crop_from_mask, crop_image, uncrop_image
from pymacruise.topology import (
    TopologyCorrectionResult,
    topology_correct,
    topology_numbers_618,
)
from pymacruise.volume_ops import VolumeGeometryError, calculate_volume
from pymacruise.workflow import WorkflowConfig, create_initial_manifest, write_run_manifest
from pymacruise.n4 import N4Config, N4DependencyError, n4_bias_correct
from pymacruise.spatial_transfer import transfer_to_reference_space

__all__ = [
    "ImageGeometryError",
    "RawPreparationResult",
    "VolumeGeometryError",
    "ValidatedInputs",
    "WorkflowInputError",
    "WorkflowProfile",
    "WorkflowConfig",
    "N4Config",
    "N4DependencyError",
    "calculate_volume",
    "cortical_separation",
    "clean_csf_boundary",
    "clean_csf_boundary_image",
    "combine_multi_atlas_toads",
    "combine_multi_atlas_toads_images",
    "merge_csf_boundary",
    "merge_csf_boundary_images",
    "multi_atlas_csf_boundary",
    "multi_atlas_csf_boundary_image",
    "n4_bias_correct",
    "transfer_to_reference_space",
    "CropParameters",
    "TopologyCorrectionResult",
    "crop_from_mask",
    "crop_image",
    "uncrop_image",
    "topology_correct",
    "topology_numbers_618",
    "create_initial_manifest",
    "load_nifti_ras",
    "save_nifti_ras",
    "prepare_raw",
    "validate_workflow_files",
    "validate_workflow_inputs",
    "write_run_manifest",
]
