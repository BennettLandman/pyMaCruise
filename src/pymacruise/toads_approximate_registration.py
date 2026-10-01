"""Atlas-aware approximate-membership method and shared registration kernels.

The lesion-descriptor initializer here is validated as an isolated overload.
No selected VandyCRUISE caller was found for a lesion-atlas approximate
registration pyramid, so no such workflow composition is exposed here. The
selected default-atlas initial pyramid lives in
``toads_default_approximate_registration``. These shared coefficient,
sampling, pooling, and line-search kernels are reused there.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, place
statistical atlas registration within the TOADS framework. This module supplies
source-derived registration kernels and one isolated atlas-aware lesion
initializer used by bounded adapters; the selected-source call order and
parameter choices are not paper claims. The kernels do not establish the
papers' full registration, iterative segmentation, or a reachable
lesion-atlas approximate-registration workflow.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_registration import (
    RegistrationCoefficientSums,
    RegistrationLineSearchResult,
    rigid_image_to_template_jacobian,
    trilinear_prior_derivatives,
)
from .toads_registration_sampling import (
    ToadsRigidRegistrationIteration,
    ToadsRigidRegistrationResult,
    sample_toad_shape_prior,
    source_block_mean_subsample,
)
from .toads_lesion_atlas_type_membership import LESION_TOADS_CLASS_TYPES


@dataclass(frozen=True)
class ToadsApproximateMembershipInitialization:
    """Sparse top-k class memberships and selected class-index slots."""

    memberships: NDArray[np.float32]
    best_class_indices: NDArray[np.uint8]
    class_types: tuple[str, ...]


def _argmax_first(values: NDArray[np.float32]) -> int:
    best = 0
    for index in range(1, values.size):
        if values[index] > values[best]:
            best = index
    return best


def _top_k(values: NDArray[np.float32], count: int, *, first: int | None = None) -> tuple[NDArray[np.uint8], NDArray[np.float32]]:
    """Match selected Numerics.bestIndex strict-tie and destructive scratch rules."""
    working = np.asarray(values, dtype=np.float32).copy()
    indices = np.empty(count, dtype=np.uint8)
    selected_values = np.empty(count, dtype=np.float32)
    for slot in range(count):
        if slot == 0 and first is not None:
            winner = int(first)
            for candidate in range(working.size):
                if working[candidate] > working[winner]:
                    winner = candidate
        else:
            winner = 0
            for candidate in range(1, working.size):
                if working[candidate] > working[winner]:
                    winner = candidate
        indices[slot] = winner
        selected_values[slot] = working[winner]
        working[winner] = np.float32(-np.inf)
    return indices, selected_values


def initialize_toads_approximate_memberships(
    image_values: ArrayLike,
    centroids: ArrayLike,
    transformed_shapes: ArrayLike,
    class_types: Sequence[str],
    *,
    approx: int,
    process_mask: ArrayLike | None = None,
) -> ToadsApproximateMembershipInitialization:
    """Run selected ``initialApproxMemberships`` on a fixed RAS subject grid.

    For each processed voxel and class, cost is the channel-ordered squared
    intensity residual plus the sum of every other class prior squared. The
    normalized reciprocal scores are sparsified to ``approx`` slots. Slot 0
    is forced to the largest transformed shape prior (strict ties choose the
    lowest class index); remaining slots are selected by strict descending
    membership, again retaining the lowest class index on ties.
    """
    image = np.asarray(image_values, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    types = tuple(str(item).upper() for item in class_types)
    if image.ndim != 4 or any(size <= 0 for size in image.shape[:3]) or image.shape[3] <= 0:
        raise ValueError("image_values must have shape (x,y,z,channel) with nonempty axes")
    if shapes.ndim != 4 or shapes.shape[:3] != image.shape[:3] or shapes.shape[-1] != len(types):
        raise ValueError("transformed_shapes must match the subject grid and class order")
    if centers.shape != (image.shape[-1], len(types)):
        raise ValueError("centroids must have shape (channel,class)")
    if types != LESION_TOADS_CLASS_TYPES:
        raise ValueError("class_types must match the selected lesion 13-class descriptor")
    if not np.isfinite(image).all() or not np.isfinite(shapes).all() or not np.isfinite(centers).all():
        raise ValueError("image, atlas, and centroid values must be finite")
    if isinstance(approx, bool) or not isinstance(approx, (int, np.integer)) or not (1 < approx < len(types)):
        raise ValueError("approx must select between two and class_count-1 candidates")
    if process_mask is None:
        process = np.ones(image.shape[:3], dtype=bool)
    else:
        process = np.asarray(process_mask, dtype=bool)
        if process.shape != image.shape[:3]:
            raise ValueError("process_mask must match the subject grid")

    slot_memberships = np.zeros(image.shape[:3] + (int(approx),), dtype=np.float32)
    candidates = np.zeros(image.shape[:3] + (int(approx),), dtype=np.uint8)
    for x in range(image.shape[0]):
        for y in range(image.shape[1]):
            for z in range(image.shape[2]):
                if not process[x, y, z]:
                    continue
                costs = np.zeros(len(types), dtype=np.float32)
                for k in range(len(types)):
                    total = np.float32(0.0)
                    for c in range(image.shape[3]):
                        residual = np.float32(image[x, y, z, c] - centers[c, k])
                        total = np.float32(total + np.float32(residual * residual))
                    for other in range(len(types)):
                        if other != k:
                            value = shapes[x, y, z, other]
                            total = np.float32(total + np.float32(value * value))
                    costs[k] = total
                reciprocals = np.empty_like(costs)
                denominator = np.float32(0.0)
                for k, cost in enumerate(costs):
                    reciprocals[k] = np.float32(1.0 / cost) if cost > np.float32(1e-20) else np.float32(1e20)
                    denominator = np.float32(denominator + reciprocals[k])
                if denominator > np.float32(0.0):
                    probabilities = np.divide(reciprocals, denominator, dtype=np.float32)
                else:
                    probabilities = reciprocals
                first = _argmax_first(shapes[x, y, z])
                ids, values = _top_k(probabilities, int(approx), first=first)
                candidates[x, y, z] = ids
                slot_memberships[x, y, z] = values
    return ToadsApproximateMembershipInitialization(slot_memberships, candidates, types)


def subsample_toads_approximate_memberships(
    memberships: ArrayLike,
    best_class_indices: ArrayLike,
    *,
    classes: int,
    factor: int,
) -> ToadsApproximateMembershipInitialization:
    """Pool top-k slots by class within blocks, then select the next top-k.

    This matches ``MultipleToadDeformableAtlas.subsampleApproxMemberships``:
    each retained full-grid slot contributes membership/factor³ to its class
    total; the coarse slots are then selected by Java strict ``bestIndex``.
    Partial blocks at the high edge are discarded.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    best = np.asarray(best_class_indices)
    if mem.ndim != 4 or best.shape != mem.shape or not np.issubdtype(best.dtype, np.integer):
        raise ValueError("memberships and best_class_indices must be matching 4D slot maps")
    if isinstance(factor, bool) or not isinstance(factor, (int, np.integer)) or factor < 2:
        raise ValueError("factor must be an integer >= 2")
    if isinstance(classes, bool) or not isinstance(classes, (int, np.integer)) or classes < 3 or classes > 127:
        raise ValueError("classes must be an integer from 3 through 127")
    if np.any((best < 0) | (best >= classes)) or not np.isfinite(mem).all():
        raise ValueError("candidate indices or memberships are invalid")
    coarse_dims = tuple(size // int(factor) for size in mem.shape[:3])
    if any(size == 0 for size in coarse_dims):
        raise ValueError("factor is too large for membership dimensions")
    scale = np.float32(int(factor) ** 3)
    coarse_mem = np.zeros(coarse_dims + (mem.shape[-1],), dtype=np.float32)
    coarse_best = np.zeros(coarse_dims + (mem.shape[-1],), dtype=np.uint8)
    for x in range(coarse_dims[0]):
        for y in range(coarse_dims[1]):
            for z in range(coarse_dims[2]):
                totals = np.zeros(int(classes), dtype=np.float32)
                for dx in range(factor):
                    for dy in range(factor):
                        for dz in range(factor):
                            sx, sy, sz = x * factor + dx, y * factor + dy, z * factor + dz
                            for slot in range(mem.shape[-1]):
                                klass = int(best[sx, sy, sz, slot])
                                contribution = np.float32(mem[sx, sy, sz, slot] / scale)
                                totals[klass] = np.float32(totals[klass] + contribution)
                ids, values = _top_k(totals, mem.shape[-1])
                coarse_best[x, y, z] = ids
                coarse_mem[x, y, z] = values
    return ToadsApproximateMembershipInitialization(coarse_mem, coarse_best, tuple("" for _ in range(classes)))


def _expand_slots(state: ToadsApproximateMembershipInitialization, classes: int) -> NDArray[np.float32]:
    expanded = np.zeros(state.memberships.shape[:3] + (classes,), dtype=np.float32)
    for index in np.ndindex(state.memberships.shape[:3]):
        seen: set[int] = set()
        for slot in range(state.memberships.shape[-1]):
            klass = int(state.best_class_indices[index + (slot,)])
            if klass in seen:
                raise ValueError("top-k candidate slots must contain unique class indices")
            seen.add(klass)
            expanded[index + (klass,)] = state.memberships[index + (slot,)]
    return expanded


def _approx_energy_and_coefficients(
    mem: NDArray[np.float32], best: NDArray[np.uint8], shapes: NDArray[np.float32],
    transform: NDArray[np.float32], *, image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float], full_image_dims: tuple[int, int, int],
    full_template_dims: tuple[int, int, int], scale: int, offset: int, subsample: int,
    registered: NDArray[np.bool_],
) -> tuple[np.float32, RegistrationCoefficientSums]:
    dims = mem.shape[:3]
    image_center = tuple(np.float32(size / 2.0) for size in full_image_dims)
    template_center = tuple(np.float32(size / 2.0) for size in full_template_dims)
    priors = np.zeros(dims + (shapes.shape[-1],), dtype=np.float32)
    derivatives = np.zeros(dims + (shapes.shape[-1], 3), dtype=np.float32)
    for x in range(dims[0]):
        for y in range(dims[1]):
            for z in range(dims[2]):
                for klass in range(1, shapes.shape[-1]):
                    sample = sample_toad_shape_prior(
                        shapes[..., klass], (x, y, z), image_dimensions_xyz=dims,
                        image_resolutions_mm=image_resolutions_mm,
                        template_resolutions_mm=template_resolutions_mm,
                        transform=tuple(float(v) for v in transform),
                        class_index=klass, bbox_low_xyz=(0, 0, 0), bbox_high_xyz=dims,
                        shape_slope=1.0, scale=scale,
                        image_center_xyz=image_center, template_center_xyz=template_center,
                    )
                    priors[x, y, z, klass] = sample.value
                    if sample.template_xyz is not None:
                        derivatives[x, y, z, klass] = trilinear_prior_derivatives(
                            shapes[..., klass], sample.template_xyz
                        )

    cost = np.float32(0.0)
    norm = np.float32(0.0)
    gradient = np.zeros(6, dtype=np.float32)
    hessian = np.zeros(6, dtype=np.float32)
    norm_gradient = np.zeros(6, dtype=np.float32)
    norm_hessian = np.zeros(6, dtype=np.float32)
    for x in range(offset, dims[0], subsample):
        for y in range(offset, dims[1], subsample):
            for z in range(offset, dims[2], subsample):
                jac = rigid_image_to_template_jacobian(
                    (x, y, z), transform, image_center,
                    image_resolutions_mm, template_resolutions_mm, scale=scale,
                )
                for slot in range(mem.shape[-1]):
                    klass = int(best[x, y, z, slot])
                    if klass == 0 or not registered[klass]:
                        continue
                    prior = priors[x, y, z, klass]
                    if prior <= np.float32(0.0):
                        continue
                    membership = mem[x, y, z, slot]
                    membership2 = np.float32(membership * membership)
                    prior2 = np.float32(prior * prior)
                    cost = np.float32(cost + np.float32(membership2 * prior2))
                    norm = np.float32(norm + prior2)
                    dp = np.zeros(6, dtype=np.float32)
                    for parameter in range(6):
                        value = np.float32(derivatives[x, y, z, klass, 0] * jac[0, parameter])
                        value = np.float32(value + np.float32(derivatives[x, y, z, klass, 1] * jac[1, parameter]))
                        dp[parameter] = np.float32(value + np.float32(derivatives[x, y, z, klass, 2] * jac[2, parameter]))
                    for parameter in range(6):
                        grad_term = np.float32(np.float32(np.float32(np.float32(2.0) * membership2) * prior) * dp[parameter])
                        hess_term = np.float32(np.float32(np.float32(np.float32(2.0) * membership2) * dp[parameter]) * dp[parameter])
                        ngrad_term = np.float32(np.float32(np.float32(2.0) * prior) * dp[parameter])
                        nhess_term = np.float32(np.float32(np.float32(2.0) * dp[parameter]) * dp[parameter])
                        gradient[parameter] = np.float32(gradient[parameter] + grad_term)
                        hessian[parameter] = np.float32(hessian[parameter] + hess_term)
                        norm_gradient[parameter] = np.float32(norm_gradient[parameter] + ngrad_term)
                        norm_hessian[parameter] = np.float32(norm_hessian[parameter] + nhess_term)
    with np.errstate(divide="ignore", invalid="ignore"):
        energy = np.float32(cost / norm)
        hessian = np.divide(hessian, norm, dtype=np.float32)
        gradient = np.divide(gradient, norm, dtype=np.float32)
        norm_gradient = np.divide(norm_gradient, norm, dtype=np.float32)
        norm_hessian = np.divide(norm_hessian, norm, dtype=np.float32)
    return energy, RegistrationCoefficientSums(cost, norm, gradient, hessian, norm_gradient, norm_hessian)


def _approx_registration_level(
    state: ToadsApproximateMembershipInitialization,
    atlas_shapes: NDArray[np.float32],
    transform: NDArray[np.float32], *,
    full_image_dims: tuple[int, int, int], full_template_dims: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float], template_resolutions_mm: tuple[float, float, float],
    scale: int, registered: NDArray[np.bool_], max_iterations: int,
    energy_tolerance: float, minimum_lambda: float,
) -> tuple[ToadsRigidRegistrationResult, tuple, tuple]:
    mem, best = state.memberships, state.best_class_indices
    atlas_level = source_block_mean_subsample(atlas_shapes, scale)
    if atlas_level.shape[:3] != tuple(size // scale for size in full_template_dims):
        raise ValueError("coarse atlas dimensions do not match source floor truncation")

    def evaluate(parameters: NDArray[np.float32]) -> np.float32:
        energy, _ = _approx_energy_and_coefficients(
            mem, best, atlas_level, parameters,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            full_image_dims=full_image_dims, full_template_dims=full_template_dims,
            scale=scale, offset=0, subsample=2, registered=registered,
        )
        return energy

    current = transform.copy()
    lam = np.float32(1.0)
    tolerance = np.float32(energy_tolerance)
    lambda_floor = np.float32(minimum_lambda)
    difference = np.float32(1.0)
    current_energy = evaluate(current)
    if not np.isfinite(current_energy):
        raise ValueError("approximate registration energy is not finite")
    history: list[ToadsRigidRegistrationIteration] = []
    coefficient_history: list[tuple[NDArray[np.float32], NDArray[np.float32]]] = []
    trial_history: list[tuple[tuple[NDArray[np.float32], np.float32], ...]] = []
    stop_reason = "iteration_limit"
    for _ in range(max_iterations):
        initial_energy, coefficients = _approx_energy_and_coefficients(
            mem, best, atlas_level, current,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            full_image_dims=full_image_dims, full_template_dims=full_template_dims,
            scale=scale, offset=0, subsample=2, registered=registered,
        )
        if not np.isfinite(coefficients.hessian).all() or not np.isfinite(coefficients.gradient).all():
            raise ValueError("approximate registration coefficients are not finite")
        coefficient_history.append((coefficients.gradient.copy(), coefficients.hessian.copy()))
        result, trials = _approx_line_search(
            current, coefficients.gradient, coefficients.hessian, initial_energy,
            evaluate, initial_lambda=lam,
        )
        trial_history.append(trials)
        current = result.transform.copy()
        lam = result.lambda_value
        difference = result.relative_change
        current_energy = result.energy
        history.append(ToadsRigidRegistrationIteration(
            current.copy(), current_energy, lam, difference,
            result.accepted, result.rejected_trials,
        ))
        if difference <= tolerance:
            stop_reason = "energy_tolerance"
            break
        if lam <= lambda_floor:
            stop_reason = "lambda_floor"
            break
    return (
        ToadsRigidRegistrationResult(current.copy(), current_energy, tuple(history), stop_reason),
        tuple(coefficient_history), tuple(trial_history),
    )


def _approx_line_search(
    transform: NDArray[np.float32], gradient: NDArray[np.float32], hessian: NDArray[np.float32],
    initial_energy: np.float32, evaluate, *, initial_lambda: np.float32,
) -> tuple[RegistrationLineSearchResult, tuple[tuple[NDArray[np.float32], np.float32], ...]]:
    current = transform.copy()
    lam = np.float32(initial_lambda)
    energy0 = np.float32(initial_energy)
    rejected = 0
    accepted = False
    energy = energy0
    trials: list[tuple[NDArray[np.float32], np.float32]] = []
    while True:
        trial = np.empty_like(current)
        for index in range(6):
            denominator = np.float32(max(np.float32(1e-30), hessian[index]))
            step = np.float32(lam / denominator)
            step = np.float32(step * gradient[index])
            trial[index] = np.float32(current[index] + step)
        energy = np.float32(evaluate(trial))
        trials.append((trial.copy(), energy))
        if energy > energy0:
            current = trial
            lam = np.float32(lam * np.float32(1.5))
            accepted = True
            break
        lam = np.float32(lam / np.float32(1.5))
        rejected += 1
        if rejected > 10:
            break
    with np.errstate(divide="ignore", invalid="ignore"):
        relative = np.float32(np.float32(energy - energy0) / energy)
    return RegistrationLineSearchResult(current, energy, lam, relative, accepted, rejected), tuple(trials)
