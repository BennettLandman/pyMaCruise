"""Later atlas-aware top-k membership to single-scale TOADS registration.

This covers the selected ``AlgorithmMultipleToads`` handoff after a
segmentation state exists: recompute atlas-aware candidates, retain the
current rigid transform through ``updateShapeRegistration``, and optionally
run ``registerShapes`` at scale one. The membership initializer and shared
sampling/coefficient kernel are reused; this is not a segmentation iteration.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, place
atlas registration within the TOADS segmentation framework. This module
implements one selected VandyCRUISE later-registration handoff after
segmentation state exists; its chosen atlas-aware overload, scale-one
optimizer route, and transform carry-forward are source-specific. It is not
the papers' complete registration or segmentation procedure.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_approximate_registration import _approx_energy_and_coefficients
from .toads_default_approximate_membership import (
    ToadsDefaultApproximateMembershipInitialization,
    initialize_toads_default_approximate_memberships,
)
from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from .toads_registration_sampling import ToadsRigidRegistrationIteration, ToadsRigidRegistrationResult
from .toads_shapes import ToadShapeBounds, _rotation_matrix, compute_toad_shape_bounds


@dataclass(frozen=True)
class ToadsDefaultSingleScaleRegistration:
    """Atlas-aware memberships and the selected optional scale-one update."""

    initialization: ToadsDefaultApproximateMembershipInitialization
    registration: ToadsRigidRegistrationResult | None
    transform: NDArray[np.float32]
    rotation: NDArray[np.float32]
    shape_bounds: tuple[ToadShapeBounds, ...]
    call_order: tuple[str, ...]
    iteration_coefficients: tuple[tuple[NDArray[np.float32], NDArray[np.float32]], ...]
    iteration_trials: tuple[tuple[tuple[NDArray[np.float32], np.float32], ...], ...]


def initialize_and_register_toads_default_single_scale(
    image_values: Any,
    centroids: Any,
    transformed_shapes: Any,
    atlas_shapes: Any,
    *,
    approx: int,
    initial_transform: Any,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    max_iterations: int = 1,
    energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
    align_shapes: bool = True,
) -> ToadsDefaultSingleScaleRegistration:
    """Compose the selected atlas-aware initialization and single-scale call.

    Inputs are already preprocessed/cropped arrays in RAS voxel order. The
    transformed subject-grid atlas maps correspond to the incoming transform;
    ``atlas_shapes`` remains on the template grid. Zero/nonpositive Hessian
    coordinates are rejected because this finite synthetic composition
    compares the source's direct ``lambda / hessian`` update without replacing
    its IEEE behavior with a denominator floor.
    """
    image = np.asarray(image_values, dtype=np.float32)
    shapes_on_subject = np.asarray(transformed_shapes, dtype=np.float32)
    shapes_on_template = np.asarray(atlas_shapes, dtype=np.float32)
    transform = np.asarray(initial_transform, dtype=np.float32)
    if image.ndim != 4 or image.shape[-1] < 1:
        raise ValueError("image_values must have shape (x,y,z,channel)")
    if shapes_on_subject.shape != image.shape[:3] + (11,):
        raise ValueError("transformed_shapes must match the subject grid and 11 classes")
    if shapes_on_template.ndim != 4 or shapes_on_template.shape[-1] != 11:
        raise ValueError("atlas_shapes must be a template grid with eleven classes")
    if not np.isfinite(image).all() or not np.isfinite(shapes_on_subject).all() or not np.isfinite(shapes_on_template).all():
        raise ValueError("image and shape inputs must be finite")
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")
    for name, resolution in (("image_resolutions_mm", image_resolutions_mm), ("template_resolutions_mm", template_resolutions_mm)):
        values = np.asarray(resolution, dtype=np.float32)
        if values.shape != (3,) or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError(f"{name} must contain three positive finite RAS voxel sizes")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    if not np.isfinite(energy_tolerance) or energy_tolerance < 0:
        raise ValueError("energy_tolerance must be finite and nonnegative")
    if not np.isfinite(minimum_lambda) or minimum_lambda < 0:
        raise ValueError("minimum_lambda must be finite and nonnegative")

    initialized = initialize_toads_default_approximate_memberships(
        image, centroids, shapes_on_subject, approx=approx,
    )
    call_order = (
        "preprocess.initialMemberships(memberships,best,atlas,approx,centroids,classes)",
        'atlas.updateShapeRegistration(memberships,best,segmentation,modality,"rigid",alignIter,1,...)',
    )
    current = transform.copy()
    registered = np.asarray([kind != "MASK" for kind in DEFAULT_ATLAS_CLASS_TYPES], dtype=bool)
    bounds = tuple(
        compute_toad_shape_bounds(
            shapes_on_template[..., klass],
            image_dimensions_xyz=image.shape[:3],
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(v) for v in current),
        )
        for klass in range(11)
    )
    if not align_shapes:
        return ToadsDefaultSingleScaleRegistration(
            initialized, None, current, _rotation_matrix(current[:3]), bounds, call_order, (), (),
        )

    call_order += ("atlas.registerShapes() -> registerAllShapes()",)
    current_energy, _ = _approx_energy_and_coefficients(
        initialized.memberships, initialized.best_class_indices, shapes_on_template, current,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        full_image_dims=image.shape[:3], full_template_dims=shapes_on_template.shape[:3],
        scale=1, offset=0, subsample=2, registered=registered,
    )
    if not np.isfinite(current_energy):
        raise ValueError("single-scale registration energy is not finite")
    lam = np.float32(1.0)
    diff = np.float32(1.0)
    history: list[ToadsRigidRegistrationIteration] = []
    coefficient_history: list[tuple[NDArray[np.float32], NDArray[np.float32]]] = []
    trial_history: list[tuple[tuple[NDArray[np.float32], np.float32], ...]] = []
    stop_reason = "iteration_limit"

    for _ in range(int(max_iterations)):
        initial_energy, coefficients = _approx_energy_and_coefficients(
            initialized.memberships, initialized.best_class_indices, shapes_on_template, current,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            full_image_dims=image.shape[:3], full_template_dims=shapes_on_template.shape[:3],
            scale=1, offset=0, subsample=2, registered=registered,
        )
        if not np.isfinite(coefficients.gradient).all() or not np.isfinite(coefficients.hessian).all():
            raise ValueError("single-scale registration coefficients are not finite")
        if np.any(coefficients.hessian <= 0):
            raise ValueError("single-scale source-faithful step requires positive Hessian coordinates")
        coefficient_history.append((coefficients.gradient.copy(), coefficients.hessian.copy()))

        trial_trace: list[tuple[NDArray[np.float32], np.float32]] = []
        trial_energy = initial_energy
        accepted = False
        rejected = 0
        candidate = current.copy()
        while True:
            trial = np.empty_like(current)
            for index in range(6):
                step = np.float32(np.float32(lam / coefficients.hessian[index]) * coefficients.gradient[index])
                trial[index] = np.float32(current[index] + step)
            trial_energy, _ = _approx_energy_and_coefficients(
                initialized.memberships, initialized.best_class_indices, shapes_on_template, trial,
                image_resolutions_mm=image_resolutions_mm,
                template_resolutions_mm=template_resolutions_mm,
                full_image_dims=image.shape[:3], full_template_dims=shapes_on_template.shape[:3],
                scale=1, offset=0, subsample=2, registered=registered,
            )
            trial_trace.append((trial.copy(), np.float32(trial_energy)))
            if trial_energy > initial_energy:
                candidate = trial
                lam = np.float32(lam * np.float32(1.5))
                accepted = True
                break
            lam = np.float32(lam / np.float32(1.5))
            rejected += 1
            if rejected > 10:
                break
        trial_history.append(tuple(trial_trace))

        with np.errstate(divide="ignore", invalid="ignore"):
            diff = np.float32(np.float32(trial_energy - initial_energy) / trial_energy)
        current = candidate
        current_energy = np.float32(trial_energy)
        history.append(ToadsRigidRegistrationIteration(
            current.copy(), current_energy, lam, diff, accepted, rejected,
        ))
        if diff <= np.float32(energy_tolerance):
            stop_reason = "energy_tolerance"
            break
        if lam <= np.float32(minimum_lambda):
            stop_reason = "lambda_floor"
            break

    result = ToadsRigidRegistrationResult(
        current.copy(), current_energy, tuple(history), stop_reason,
    )
    final_bounds = tuple(
        compute_toad_shape_bounds(
            shapes_on_template[..., klass],
            image_dimensions_xyz=image.shape[:3],
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(v) for v in current),
        )
        for klass in range(11)
    )
    return ToadsDefaultSingleScaleRegistration(
        initialized, result, current.copy(), _rotation_matrix(current[:3]), final_bounds,
        call_order, tuple(coefficient_history), tuple(trial_history),
    )
