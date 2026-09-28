"""Source-ordered operations for TOADS rigid atlas registration."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def exact_membership_registration_energy(
    memberships: ArrayLike,
    shape_priors: ArrayLike,
    *,
    registered_classes: ArrayLike | None = None,
    offset: int = 0,
    subsample: int = 1,
) -> np.float32:
    """Evaluate the selected TOADS exact-membership rigid-registration energy.

    This is the source's `computeExactRegistrationEnergy` accumulation after
    priors have been sampled. The first class (background) is always skipped.
    Only positive priors contribute, and the energy is the sum of squared
    membership-prior products divided by the sum of squared priors. The
    source uses float32 accumulation in x/y/z/class loop order.

    Transform mapping, closest-edge interpolation, and prior derivatives are
    separate operations. This function is an objective kernel, not the atlas
    optimizer or a complete registration substitute.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    prior = np.asarray(shape_priors, dtype=np.float32)
    if mem.ndim != 4 or mem.shape != prior.shape:
        raise ValueError("memberships and shape_priors must share a 4D (x, y, z, class) grid")
    if any(size < 1 for size in mem.shape[:3]) or mem.shape[3] < 2:
        raise ValueError("registration arrays need nonempty spatial axes and at least two classes")
    if not np.isfinite(mem).all() or not np.isfinite(prior).all():
        raise ValueError("memberships and shape_priors must be finite")
    if not isinstance(offset, (int, np.integer)) or offset < 0:
        raise ValueError("offset must be a non-negative integer")
    if not isinstance(subsample, (int, np.integer)) or subsample < 1:
        raise ValueError("subsample must be a positive integer")
    if offset >= min(mem.shape[:3]):
        raise ValueError("offset must be inside every spatial axis")

    if registered_classes is None:
        registered = np.ones(mem.shape[3], dtype=bool)
    else:
        registered = np.asarray(registered_classes, dtype=bool)
        if registered.shape != (mem.shape[3],):
            raise ValueError("registered_classes must have one value per class")

    cost = np.float32(0.0)
    norm = np.float32(0.0)
    for x in range(offset, mem.shape[0], subsample):
        for y in range(offset, mem.shape[1], subsample):
            for z in range(offset, mem.shape[2], subsample):
                for k in range(1, mem.shape[3]):
                    if not registered[k]:
                        continue
                    prior_value = prior[x, y, z, k]
                    if prior_value > np.float32(0.0):
                        membership = mem[x, y, z, k]
                        membership_square = np.float32(membership * membership)
                        prior_square = np.float32(prior_value * prior_value)
                        contribution = np.float32(membership_square * prior_square)
                        cost = np.float32(cost + contribution)
                        norm = np.float32(norm + prior_square)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.float32(cost / norm)


def exact_registration_measure_coefficients(
    gradient: ArrayLike,
    hessian: ArrayLike,
    norm_gradient: ArrayLike,
    norm_hessian: ArrayLike,
    norm: float,
    *,
    source_faithful: bool = True,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Apply the selected source's final rigid-registration coefficient step.

    The VandyCRUISE method first assigns ``hessian[i] = hessian[i] / norm``
    and then calls ``registrationMeasureGradient(hessian[i], nhess[i],
    norm, normalize)``. The helper returns its first argument divided by
    ``norm``, so this source call uses the normalized Hessian as the gradient
    numerator and divides by ``norm`` again. ``source_faithful=True`` preserves
    this observed behavior. Set it false to compute the apparent intended
    ``gradient / norm`` for experiments; that corrected behavior is not the
    historical default and has no output-equivalence claim.

    The source's current helper implementations ignore ``norm_gradient``,
    ``norm_hessian``, and the ``normalize`` scaling argument.
    """
    grad = np.asarray(gradient, dtype=np.float32)
    hess = np.asarray(hessian, dtype=np.float32)
    ngrad = np.asarray(norm_gradient, dtype=np.float32)
    nhess = np.asarray(norm_hessian, dtype=np.float32)
    if grad.ndim != 1 or hess.shape != grad.shape:
        raise ValueError("gradient and hessian must be matching one-dimensional arrays")
    if ngrad.shape != grad.shape or nhess.shape != grad.shape:
        raise ValueError("norm gradient and Hessian must match the coefficient arrays")
    if not all(np.isfinite(values).all() for values in (grad, hess, ngrad, nhess)):
        raise ValueError("registration coefficients must be finite")
    norm_value = np.float32(norm)
    if not np.isfinite(norm_value) or norm_value < 0:
        raise ValueError("norm must be finite and non-negative")
    with np.errstate(divide="ignore", invalid="ignore"):
        scaled_hessian = np.divide(hess, norm_value, dtype=np.float32)
        if source_faithful:
            scaled_gradient = np.divide(scaled_hessian, norm_value, dtype=np.float32)
        else:
            scaled_gradient = np.divide(grad, norm_value, dtype=np.float32)
    return scaled_gradient, scaled_hessian


def synthetic_registration_objective_inputs() -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Return the deterministic 5×4×3×4 Java-reference objective phantom."""
    x, y, z, k = np.indices((5, 4, 3, 4), dtype=np.int32)
    memberships = (((7 * x + 3 * y + 5 * z + 2 * k) % 19) / np.float32(18.0)).astype(np.float32)
    priors = (((2 * x + 5 * y + 3 * z + 4 * k) % 11) / np.float32(10.0)).astype(np.float32)
    return memberships, priors
