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

from pymacruise.atlas_em_preprocessing import (
    AtlasEMCroppedImages,
    AtlasEMNormalizedImages,
    recrop_atlas_em_after_registration,
)
from pymacruise.toads_registration_sampling import (
    ToadsRigidRegistrationPyramidResult,
    ToadsRegistrationPyramidLevel,
    build_toads_registration_pyramid,
    register_toads_atlas_rigid_pyramid,
)


@dataclass(frozen=True)
class AtlasEMRegistrationRecropResult:
    """Rigid-registration history and the source-ordered following crop."""

    registration: ToadsRigidRegistrationPyramidResult
    recropped: AtlasEMCroppedImages


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
