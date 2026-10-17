"""AtlasEM rigid-registration preparation backed by shared TOADS kernels.

AtlasEM's ``DemonToadDeformableAtlas.registerShapesPyramid`` builds its
membership and per-class shape pyramids with ``ImageFunctions.subsample``.
The array operation is shared with the already source-checked TOADS pyramid
builder; the optimizer and transform update remain separate AtlasEM work.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pymacruise.atlas_em_descriptor import (
    AtlasEMRegistrationAtlasInputs,
    load_atlas_em_registration_inputs,
    parse_atlas_em_descriptor,
)
from pymacruise.atlas_em_preprocessing import (
    AtlasEMDescriptorInitialization,
    AtlasEMCroppedImages,
    AtlasEMNormalizedImages,
    initialize_atlas_em_memberships_from_descriptor,
    recrop_atlas_em_after_registration,
)
from pymacruise.toads_registration_sampling import (
    ToadsRigidRegistrationResult,
    ToadsRigidRegistrationPyramidResult,
    ToadsRegistrationPyramidLevel,
    build_toads_registration_pyramid,
    register_toads_atlas_rigid,
    register_toads_atlas_rigid_pyramid,
)

ATLAS_EM_PLUGIN_DEFAULT_LEVELS = 4
ATLAS_EM_PLUGIN_DEFAULT_INITIAL_REGISTRATION_ITERATIONS = 50
ATLAS_EM_PLUGIN_DEFAULT_POST_RECROP_REGISTRATION_ITERATIONS = 1
ATLAS_EM_PLUGIN_DEFAULT_DEMONS_SMOOTHING = 1.0
ATLAS_EM_PLUGIN_DEFAULT_DEMONS_SCALE = 2.0


@dataclass(frozen=True)
class AtlasEMRegistrationRecropResult:
    """Rigid-registration history and the source-ordered following crop."""

    registration: ToadsRigidRegistrationPyramidResult
    recropped: AtlasEMCroppedImages


@dataclass(frozen=True)
class AtlasEMDescriptorRegistrationResult:
    """Descriptor-loaded atlas and the rigid registration computed from it."""

    atlas: AtlasEMRegistrationAtlasInputs
    registration: ToadsRigidRegistrationPyramidResult


@dataclass(frozen=True)
class AtlasEMDescriptorRegistrationRecropResult:
    """Descriptor-loaded AtlasEM registration and its native-grid recrop."""

    atlas: AtlasEMRegistrationAtlasInputs
    registration: ToadsRigidRegistrationPyramidResult
    recropped: AtlasEMCroppedImages


@dataclass(frozen=True)
class AtlasEMInitialRegistrationResult:
    """Descriptor-driven memberships followed by registration and recropping."""

    initialization: AtlasEMDescriptorInitialization
    atlas: AtlasEMRegistrationAtlasInputs
    registration: ToadsRigidRegistrationPyramidResult
    recropped: AtlasEMCroppedImages


@dataclass(frozen=True)
class AtlasEMTwoStageRegistrationResult:
    """Initial pyramid registration, recrop, and post-recrop rigid update."""

    initialization: AtlasEMDescriptorInitialization
    atlas: AtlasEMRegistrationAtlasInputs
    initial_registration: ToadsRigidRegistrationPyramidResult
    recropped: AtlasEMCroppedImages
    post_recrop_initialization: AtlasEMDescriptorInitialization
    post_recrop_registration: ToadsRigidRegistrationResult


def register_atlas_em_rigid_single_scale(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: tuple[float, float, float, float, float, float],
    registration_iterations: int,
    registered_classes: ArrayLike | None = None,
) -> ToadsRigidRegistrationResult:
    """Run AtlasEM's post-recrop, full-resolution rigid registration stage.

    The selected ``AlgorithmAtlasEMSegmentation`` calls ``updateShapeRegistration``
    with one level after initializing priors and recomputing memberships, then
    calls ``registerShapes()``. ``DemonToadDeformableAtlas`` leaves voxel
    subsampling at 2 and its exact single-scale coefficient routine assigns the
    normalized Hessian into the gradient slot before normalizing again. The
    source-faithful shared registration kernel preserves that call order.

    This adapter covers the registration stage only. Membership arrays and
    transform must already reflect the preceding AtlasEM initialization and
    recrop. The Java coefficient loop samples the complete current image grid;
    its per-class shape bounds are recomputed after registration and do not
    constrain optimizer samples. Iteration count is explicit because the
    MATLAB run wrapper does not expose the Java plugin's main-registration
    default.
    """
    if isinstance(registration_iterations, bool) or not isinstance(
        registration_iterations, (int, np.integer)
    ) or registration_iterations < 1:
        raise ValueError("registration_iterations must be a positive integer")
    membership_values = np.asarray(memberships)
    shape_values = np.asarray(atlas_shapes)
    if membership_values.ndim != 4 or shape_values.ndim != 4:
        raise ValueError("memberships and atlas_shapes must be 4D arrays")
    if membership_values.shape[-1] != shape_values.shape[-1]:
        raise ValueError("memberships and atlas_shapes must contain the same classes")
    # DemonToadDeformableAtlas.computeRegistrationCoefficients() has no class
    # bounding-box check; it samples every stride-two image-grid point. Bounds
    # retained from the shared TOADS registration API must therefore span the
    # image for this source-specific adapter.
    full_low = np.zeros((membership_values.shape[-1], 3), dtype=np.int64)
    full_high = np.broadcast_to(
        np.asarray(membership_values.shape[:3], dtype=np.int64), full_low.shape
    ).copy()
    return register_toads_atlas_rigid(
        memberships,
        atlas_shapes,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=initial_transform,
        class_bounds_low_xyz=full_low,
        class_bounds_high_xyz=full_high,
        registered_classes=registered_classes,
        offset=0,
        subsample=2,
        max_iterations=int(registration_iterations),
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
        source_faithful=True,
        allow_nonfinite_trials=True,
        scale=1.0,
        # DemonToadDeformableAtlas enables matrix-first coordinates for linear
        # rigid transforms; keep its float32 rounding path in this second stage.
        precomputed_affine=True,
    )


def build_atlas_em_rigid_registration_pyramid(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    levels: int,
) -> tuple[ToadsRegistrationPyramidLevel, ...]:
    """Build AtlasEM's source-ordered coarse membership/shape levels.

    For ``levels=L``, the selected Java wrapper visits scales
    ``2**(L-1), ..., 2`` from coarse to fine and does not run a scale-1
    registration stage. Inputs are native-grid ``(x,y,z,class)`` float32
    arrays; each level is block-mean subsampled directly from those arrays.
    This only prepares registration inputs. It does not estimate a rigid
    transform or claim that the TOADS optimizer equals AtlasEM registration.
    """
    return build_toads_registration_pyramid(memberships, atlas_shapes, levels=levels)


def register_atlas_em_rigid_pyramid(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    levels: int,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: tuple[float, float, float, float, float, float],
    registration_iterations: int,
    registered_classes: ArrayLike | None = None,
) -> ToadsRigidRegistrationPyramidResult:
    """Run the selected AtlasEM rigid-registration schedule on one input.

    This connects AtlasEM's ``"rigid"``/pyramid configuration to the shared
    transform sampler and coefficient kernels. AtlasEM source configures
    offset zero, stride two, per-level lambda 1, line-search factor 1.5,
    rejection limit 10, relative-change threshold 1e-6, and minimum lambda
    0.001. Its ``DemonToadDeformableAtlas`` gradient uses the direct gradient
    divided by the prior norm, so ``source_faithful=False`` selects that
    source class's formula (the adjacent MultipleToad source has a distinct
    historical Hessian/gradient mix-up).

    The wrapper-level settings and reused numerical kernels are source-backed;
    this API is not yet validated end-to-end against a Java execution of
    ``DemonToadDeformableAtlas``. In particular, atlas descriptor loading and
    registration output parity still require that direct class comparison.
    """
    if not isinstance(registration_iterations, int) or registration_iterations < 0:
        raise ValueError("registration_iterations must be a non-negative integer")
    return register_toads_atlas_rigid_pyramid(
        memberships,
        atlas_shapes,
        levels=levels,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=initial_transform,
        registered_classes=registered_classes,
        max_iterations=registration_iterations,
        subsample_factor=2,
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
        source_faithful=False,
        allow_nonfinite_trials=True,
    )


def register_atlas_em_descriptor_pyramid(
    memberships: ArrayLike,
    descriptor_path: str,
    *,
    levels: int,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float] | None = None,
    initial_transform: tuple[float, float, float, float, float, float] = (0.0,) * 6,
    registration_iterations: int,
) -> AtlasEMDescriptorRegistrationResult:
    """Load a source-ordered `.nii.gz` atlas and run its rigid registration.

    The descriptor supplies the class order, registered-class mask, and atlas
    resolution; subject spacing/initial transform and the iteration cap remain
    explicit run inputs. This is the selected adapter composition and still
    requires same-input Java class comparison for optimizer equivalence.
    """
    descriptor = parse_atlas_em_descriptor(descriptor_path)
    atlas = load_atlas_em_registration_inputs(descriptor)
    template_resolution = (
        descriptor.shape_resolutions_mm
        if template_resolutions_mm is None else template_resolutions_mm
    )
    registration = register_atlas_em_rigid_pyramid(
        memberships,
        atlas.shape_priors,
        levels=levels,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolution,
        initial_transform=initial_transform,
        registration_iterations=registration_iterations,
        registered_classes=atlas.registered_classes,
    )
    return AtlasEMDescriptorRegistrationResult(atlas, registration)


def register_and_recrop_atlas_em(
    crop: AtlasEMCroppedImages,
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    normalized: AtlasEMNormalizedImages,
    *,
    background_ratio: float,
    levels: int,
    template_resolutions_mm: tuple[float, float, float],
    registration_iterations: int,
    registered_classes: ArrayLike | None = None,
) -> AtlasEMRegistrationRecropResult:
    """Run rigid registration, then AtlasEM's native-grid recrop handoff.

    This composes the source-ordered pyramid registration adapter with the
    Java-checked uncrop, recompute-bounds, and crop sequence. ``memberships``
    must be on the current cropped image grid; ``normalized`` retains the
    full native-grid image channels and robust bounds. The optimizer result
    is passed directly to the recrop transform update. Direct execution of
    this composed handoff against the selected Java class remains pending.
    """
    image_values = np.asarray(crop.images, dtype=np.float32)
    mem = np.asarray(memberships, dtype=np.float32)
    if image_values.ndim != 4 or mem.ndim != 4 or tuple(mem.shape[:3]) != tuple(image_values.shape[:3]):
        raise ValueError("memberships and the current crop must share a 3D grid")
    if tuple(normalized.images.shape[:3]) != crop.original_dimensions_xyz:
        raise ValueError("normalized native images must match the crop's original dimensions")
    if tuple(normalized.images.shape[3:]) != tuple(image_values.shape[3:]):
        raise ValueError("normalized and cropped image channel counts must match")
    registration = register_atlas_em_rigid_pyramid(
        mem,
        atlas_shapes,
        levels=levels,
        image_resolutions_mm=crop.signed_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=tuple(float(value) for value in crop.transform),
        registration_iterations=registration_iterations,
        registered_classes=registered_classes,
    )
    recropped = recrop_atlas_em_after_registration(
        crop,
        registration.transform,
        normalized.robust_minima,
        normalized.robust_maxima,
        background_ratio=background_ratio,
    )
    return AtlasEMRegistrationRecropResult(registration, recropped)


def register_and_recrop_atlas_em_from_descriptor(
    crop: AtlasEMCroppedImages,
    memberships: ArrayLike,
    normalized: AtlasEMNormalizedImages,
    descriptor_path: str,
    *,
    background_ratio: float,
    levels: int,
    registration_iterations: int,
    template_resolutions_mm: tuple[float, float, float] | None = None,
) -> AtlasEMDescriptorRegistrationRecropResult:
    """Load the AtlasEM atlas, register it, then perform native-grid recropping."""
    descriptor = parse_atlas_em_descriptor(descriptor_path)
    atlas = load_atlas_em_registration_inputs(descriptor)
    template_resolution = (
        descriptor.shape_resolutions_mm
        if template_resolutions_mm is None else template_resolutions_mm
    )
    handoff = register_and_recrop_atlas_em(
        crop,
        memberships,
        atlas.shape_priors,
        normalized,
        background_ratio=background_ratio,
        levels=levels,
        template_resolutions_mm=template_resolution,
        registration_iterations=registration_iterations,
        registered_classes=atlas.registered_classes,
    )
    return AtlasEMDescriptorRegistrationRecropResult(
        atlas, handoff.registration, handoff.recropped
    )


def initialize_register_recrop_atlas_em_from_descriptor(
    crop: AtlasEMCroppedImages,
    normalized: AtlasEMNormalizedImages,
    descriptor_path: str,
    modalities: tuple[str, ...],
    *,
    background_ratio: float,
    levels: int,
    registration_iterations: int,
    template_resolutions_mm: tuple[float, float, float] | None = None,
) -> AtlasEMInitialRegistrationResult:
    """Compose descriptor-driven initial memberships, registration, and recrop."""
    crop_normalized = AtlasEMNormalizedImages(
        np.asarray(crop.images, dtype=np.float32),
        normalized.robust_minima,
        normalized.robust_maxima,
    )
    initialization = initialize_atlas_em_memberships_from_descriptor(
        crop_normalized, descriptor_path, modalities
    )
    handoff = register_and_recrop_atlas_em_from_descriptor(
        crop,
        initialization.memberships,
        normalized,
        descriptor_path,
        background_ratio=background_ratio,
        levels=levels,
        registration_iterations=registration_iterations,
        template_resolutions_mm=template_resolutions_mm,
    )
    return AtlasEMInitialRegistrationResult(
        initialization,
        handoff.atlas,
        handoff.registration,
        handoff.recropped,
    )


def initialize_register_recrop_register_atlas_em_from_descriptor(
    crop: AtlasEMCroppedImages,
    normalized: AtlasEMNormalizedImages,
    descriptor_path: str,
    modalities: tuple[str, ...],
    *,
    background_ratio: float,
    levels: int = ATLAS_EM_PLUGIN_DEFAULT_LEVELS,
    initial_registration_iterations: int = ATLAS_EM_PLUGIN_DEFAULT_INITIAL_REGISTRATION_ITERATIONS,
    post_recrop_registration_iterations: int = ATLAS_EM_PLUGIN_DEFAULT_POST_RECROP_REGISTRATION_ITERATIONS,
    template_resolutions_mm: tuple[float, float, float] | None = None,
) -> AtlasEMTwoStageRegistrationResult:
    """Compose both source-ordered AtlasEM rigid-registration stages.

    The first registration runs on the initial crop and is followed by native
    grid uncrop/rebound/recrop. AtlasEM then rebuilds initial memberships on
    that recropped grid before its one-level full-resolution rigid update.
    Defaults are taken from the selected MIPAV plugin constructor call; callers
    can override both caps for reduced synthetic cases or matched settings.
    """
    first = initialize_register_recrop_atlas_em_from_descriptor(
        crop,
        normalized,
        descriptor_path,
        modalities,
        background_ratio=background_ratio,
        levels=levels,
        registration_iterations=initial_registration_iterations,
        template_resolutions_mm=template_resolutions_mm,
    )
    template_resolution = (
        parse_atlas_em_descriptor(descriptor_path).shape_resolutions_mm
        if template_resolutions_mm is None else template_resolutions_mm
    )
    recropped_normalized = AtlasEMNormalizedImages(
        np.asarray(first.recropped.images, dtype=np.float32),
        normalized.robust_minima,
        normalized.robust_maxima,
    )
    second_initialization = initialize_atlas_em_memberships_from_descriptor(
        recropped_normalized, descriptor_path, modalities
    )
    final_registration = register_atlas_em_rigid_single_scale(
        second_initialization.memberships,
        first.atlas.shape_priors,
        image_resolutions_mm=first.recropped.signed_resolutions_mm,
        template_resolutions_mm=template_resolution,
        initial_transform=tuple(float(value) for value in first.recropped.transform),
        registration_iterations=post_recrop_registration_iterations,
        registered_classes=first.atlas.registered_classes,
    )
    return AtlasEMTwoStageRegistrationResult(
        first.initialization,
        first.atlas,
        first.registration,
        first.recropped,
        second_initialization,
        final_registration,
    )
