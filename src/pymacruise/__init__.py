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
from pymacruise.ace import ACEArrayResult, ACEImageResult, ace_enhance, ace_enhance_arrays
from pymacruise.skeleton import ThinACEImageResult, thin_ace, thin_ace_image
from pymacruise.gvf import GVFImageResult, generalized_gvf, generalized_gvf_image
from pymacruise.surfaces import (
    SurfaceMesh,
    convert_surface_to_ras,
    extract_levelset_surface,
    read_legacy_vtk,
    surface_vertex_normals,
    transform_surface_ras,
    write_legacy_vtk,
)
from pymacruise.thickness import (
    SurfaceThicknessResult,
    ThicknessVolumeResult,
    closest_point_thickness,
    levelset_subtraction_thickness,
)
from pymacruise.segmentation_refine import (
    SegmentationRefinementResult,
    refine_segmentation_by_surfaces,
)
from pymacruise.registration import (
    ANTsRegistrationError,
    RegistrationConfig,
    RegistrationResult,
    apply_image_transform,
    apply_transforms_to_points_ras,
    apply_transforms_to_surface_ras,
    register_flair_to_t1,
    register_to_reference,
)
from pymacruise.toads_registration import (
    RegistrationCoefficientSums,
    RegistrationLineSearchResult,
    exact_registration_coefficient_contributions,
    exact_membership_registration_energy,
    exact_membership_registration_coefficients,
    exact_registration_gradient_descent_step,
    exact_registration_measure_coefficients,
)
from pymacruise.lesion_output_io import (
    lesion_toads_output_images,
    save_lesion_toads_outputs,
)
from pymacruise.toads_output_io import (
    save_toads_dura_wrapper_outputs,
    save_toads_outputs,
    toads_output_images,
)
from pymacruise.toads_alignment_startup import (
    ToadsAtlasAlignmentRun,
    ToadsAlignmentRegistration,
    ToadsAlignmentStartup,
    ToadsRigidRegistrationSetup,
    prepare_toads_alignment_startup,
    register_toads_alignment_startup,
    register_toads_from_atlas_assets,
)
from pymacruise.toads_atlas import (
    ToadsAtlas,
    ToadsAtlasData,
    ToadsAtlasRASData,
    load_toads_atlas,
    load_toads_atlas_data,
    load_toads_atlas_ras_data,
    load_toads_intensity_priors,
    load_toads_lesion_priors,
    toads_atlas_shapes_in_structure_order,
)
from pymacruise.toads_post_refresh_registration import (
    ToadsPostRefreshRegistration,
    register_toads_after_template_refresh,
)
from pymacruise.toads_segmentation_initialization import (
    ToadsSegmentationInitialization,
    initialize_toads_segmentation_state,
)
from pymacruise.toads_segmentation_constructor import (
    ToadsSegmentationConstructorState,
    initialize_toads_segmentation_constructor_state,
)

__all__ = [
    "ImageGeometryError",
    "ACEArrayResult",
    "ACEImageResult",
    "ThinACEImageResult",
    "GVFImageResult",
    "SurfaceMesh",
    "SurfaceThicknessResult",
    "ThicknessVolumeResult",
    "SegmentationRefinementResult",
    "RegistrationConfig",
    "RegistrationResult",
    "ANTsRegistrationError",
    "RawPreparationResult",
    "VolumeGeometryError",
    "ValidatedInputs",
    "WorkflowInputError",
    "WorkflowProfile",
    "WorkflowConfig",
    "N4Config",
    "N4DependencyError",
    "calculate_volume",
    "ace_enhance",
    "ace_enhance_arrays",
    "thin_ace",
    "thin_ace_image",
    "generalized_gvf",
    "generalized_gvf_image",
    "extract_levelset_surface",
    "convert_surface_to_ras",
    "read_legacy_vtk",
    "write_legacy_vtk",
    "surface_vertex_normals",
    "transform_surface_ras",
    "closest_point_thickness",
    "levelset_subtraction_thickness",
    "refine_segmentation_by_surfaces",
    "register_to_reference",
    "exact_membership_registration_energy",
    "exact_membership_registration_coefficients",
    "exact_registration_measure_coefficients",
    "RegistrationCoefficientSums",
    "RegistrationLineSearchResult",
    "exact_registration_coefficient_contributions",
    "exact_registration_gradient_descent_step",
    "register_flair_to_t1",
    "apply_image_transform",
    "apply_transforms_to_points_ras",
    "apply_transforms_to_surface_ras",
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
    "lesion_toads_output_images",
    "save_lesion_toads_outputs",
    "toads_output_images",
    "save_toads_outputs",
    "save_toads_dura_wrapper_outputs",
    "ToadsAlignmentStartup",
    "prepare_toads_alignment_startup",
    "ToadsAlignmentRegistration",
    "register_toads_alignment_startup",
    "ToadsAtlasAlignmentRun",
    "ToadsRigidRegistrationSetup",
    "register_toads_from_atlas_assets",
    "ToadsAtlas",
    "ToadsAtlasData",
    "ToadsAtlasRASData",
    "load_toads_atlas",
    "load_toads_atlas_data",
    "load_toads_atlas_ras_data",
    "load_toads_intensity_priors",
    "load_toads_lesion_priors",
    "toads_atlas_shapes_in_structure_order",
    "ToadsPostRefreshRegistration",
    "register_toads_after_template_refresh",
    "ToadsSegmentationInitialization",
    "initialize_toads_segmentation_state",
    "ToadsSegmentationConstructorState",
    "initialize_toads_segmentation_constructor_state",
]
