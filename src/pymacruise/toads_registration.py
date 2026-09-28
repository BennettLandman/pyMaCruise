"""Source-ordered operations for TOADS rigid atlas registration."""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class RegistrationCoefficientSums:
    """Float32 sums accumulated by the selected exact registration loop."""

    cost: np.float32
    norm: np.float32
    gradient: NDArray[np.float32]
    hessian: NDArray[np.float32]
    norm_gradient: NDArray[np.float32]
    norm_hessian: NDArray[np.float32]


@dataclass(frozen=True)
class RegistrationLineSearchResult:
    """One source-compatible TOADS registration gradient-descent step."""

    transform: NDArray[np.float32]
    energy: np.float32
    lambda_value: np.float32
    relative_change: np.float32
    accepted: bool
    rejected_trials: int


def exact_registration_gradient_descent_step(
    transform: ArrayLike,
    gradient: ArrayLike,
    hessian: ArrayLike,
    initial_energy: float,
    evaluate_energy: Callable[[NDArray[np.float32]], float],
    *,
    initial_lambda: float = 1.0,
    lambda_factor: float = 1.5,
    max_rejections: int = 10,
) -> RegistrationLineSearchResult:
    """Run one selected TOADS diagonal-gradient update with line search.

    Trials use the source expression `transform + lambda / hessian * gradient`
    in float32. The selected method maximizes energy, accepts only strict
    improvement, scales lambda by `lambda_factor`, and stops after the source
    rejection-count condition (`rejections > max_rejections`).
    """
    current = np.asarray(transform, dtype=np.float32).copy()
    grad = np.asarray(gradient, dtype=np.float32)
    hess = np.asarray(hessian, dtype=np.float32)
    if current.ndim != 1 or grad.shape != current.shape or hess.shape != current.shape:
        raise ValueError("transform, gradient, and hessian must be matching vectors")
    if current.size == 0 or not all(np.isfinite(v).all() for v in (current, grad, hess)):
        raise ValueError("registration vectors must be nonempty and finite")
    energy0 = np.float32(initial_energy)
    lam = np.float32(initial_lambda)
    factor = np.float32(lambda_factor)
    if not np.isfinite(energy0):
        raise ValueError("initial_energy must be finite")
    if not np.isfinite(lam) or lam <= 0 or not np.isfinite(factor) or factor <= 1:
        raise ValueError("initial_lambda must be positive and lambda_factor must exceed one")
    if not isinstance(max_rejections, (int, np.integer)) or max_rejections < 0:
        raise ValueError("max_rejections must be a non-negative integer")

    energy = energy0
    rejected = 0
    accepted = False
    while True:
        trial = np.empty_like(current)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            for index in range(current.size):
                step = np.float32(lam / hess[index])
                step = np.float32(step * grad[index])
                trial[index] = np.float32(current[index] + step)
        energy = np.float32(evaluate_energy(trial))
        if energy > energy0:
            current = trial
            lam = np.float32(lam * factor)
            accepted = True
            break
        lam = np.float32(lam / factor)
        rejected += 1
        if rejected > max_rejections:
            break
    with np.errstate(divide="ignore", invalid="ignore"):
        relative_change = np.float32(np.float32(energy - energy0) / energy)
    return RegistrationLineSearchResult(current, energy, lam, relative_change, accepted, rejected)


def exact_registration_coefficient_contributions(
    memberships: ArrayLike,
    shape_priors: ArrayLike,
    prior_gradients_xyz: ArrayLike,
    coordinate_jacobian: ArrayLike,
    *,
    registered_classes: ArrayLike | None = None,
) -> RegistrationCoefficientSums:
    """Accumulate one voxel's exact-membership registration coefficients.

    This follows the selected source's registered-class loop after a prior
    and its x/y/z derivatives have been sampled. It preserves float32 source
    operation order for cost/norm, gradient/Hessian, and norm derivatives.
    The caller performs final measure normalization separately.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    prior = np.asarray(shape_priors, dtype=np.float32)
    prior_grad = np.asarray(prior_gradients_xyz, dtype=np.float32)
    jacobian = np.asarray(coordinate_jacobian, dtype=np.float32)
    if mem.ndim != 1 or mem.size < 2 or prior.shape != mem.shape:
        raise ValueError("memberships and shape_priors must be matching class vectors")
    if prior_grad.shape != (mem.size, 3):
        raise ValueError("prior_gradients_xyz must have shape (classes, 3)")
    if jacobian.shape != (3, 6):
        raise ValueError("coordinate_jacobian must have shape (3, 6)")
    if not all(np.isfinite(v).all() for v in (mem, prior, prior_grad, jacobian)):
        raise ValueError("registration contribution inputs must be finite")
    if registered_classes is None:
        registered = np.ones(mem.size, dtype=bool)
    else:
        registered = np.asarray(registered_classes, dtype=bool)
        if registered.shape != mem.shape:
            raise ValueError("registered_classes must have one value per class")

    cost = np.float32(0.0)
    norm = np.float32(0.0)
    gradient = np.zeros(6, dtype=np.float32)
    hessian = np.zeros(6, dtype=np.float32)
    norm_gradient = np.zeros(6, dtype=np.float32)
    norm_hessian = np.zeros(6, dtype=np.float32)
    for k in range(1, mem.size):
        p = prior[k]
        if not registered[k] or p <= np.float32(0.0):
            continue
        m = mem[k]
        m2 = np.float32(m * m)
        p2 = np.float32(p * p)
        cost = np.float32(cost + np.float32(np.float32(m2 * p) * p))
        norm = np.float32(norm + p2)
        dp = np.empty(6, dtype=np.float32)
        for parameter in range(6):
            value = np.float32(prior_grad[k, 0] * jacobian[0, parameter])
            value = np.float32(value + np.float32(prior_grad[k, 1] * jacobian[1, parameter]))
            dp[parameter] = np.float32(value + np.float32(prior_grad[k, 2] * jacobian[2, parameter]))
        for parameter in range(6):
            gradient_term = np.float32(np.float32(2.0) * m)
            gradient_term = np.float32(gradient_term * m)
            gradient_term = np.float32(gradient_term * p)
            gradient_term = np.float32(gradient_term * dp[parameter])
            gradient[parameter] = np.float32(gradient[parameter] + gradient_term)
            hessian_term = np.float32(np.float32(2.0) * m)
            hessian_term = np.float32(hessian_term * m)
            hessian_term = np.float32(hessian_term * dp[parameter])
            hessian_term = np.float32(hessian_term * dp[parameter])
            hessian[parameter] = np.float32(hessian[parameter] + hessian_term)
            norm_gradient_term = np.float32(np.float32(np.float32(2.0) * p) * dp[parameter])
            norm_gradient[parameter] = np.float32(norm_gradient[parameter] + norm_gradient_term)
            norm_hessian_term = np.float32(np.float32(np.float32(2.0) * dp[parameter]) * dp[parameter])
            norm_hessian[parameter] = np.float32(norm_hessian[parameter] + norm_hessian_term)
    return RegistrationCoefficientSums(
        cost, norm, gradient, hessian, norm_gradient, norm_hessian
    )


def exact_membership_registration_coefficients(
    memberships: ArrayLike,
    shape_priors: ArrayLike,
    prior_gradients_xyz: ArrayLike,
    coordinate_jacobians: ArrayLike,
    *,
    registered_classes: ArrayLike | None = None,
    offset: int = 0,
    subsample: int = 1,
) -> RegistrationCoefficientSums:
    """Accumulate exact-membership rigid coefficients over a sampled grid.

    Inputs contain source-sampled priors, prior x/y/z derivatives, and a
    per-voxel `(3, 6)` rigid Jacobian. Loops follow the selected Java
    x/y/z/class ordering so each running sum retains float32 source semantics.
    This numerical kernel does not perform coordinate mapping or atlas
    interpolation.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    prior = np.asarray(shape_priors, dtype=np.float32)
    prior_grad = np.asarray(prior_gradients_xyz, dtype=np.float32)
    jacobian = np.asarray(coordinate_jacobians, dtype=np.float32)
    if mem.ndim != 4 or mem.shape != prior.shape or mem.shape[-1] < 2:
        raise ValueError("membership and prior arrays must share a 4D grid with at least two classes")
    if prior_grad.shape != mem.shape + (3,):
        raise ValueError("prior_gradients_xyz must have shape (x, y, z, class, 3)")
    if jacobian.shape != mem.shape[:3] + (3, 6):
        raise ValueError("coordinate_jacobians must have shape (x, y, z, 3, 6)")
    if not all(np.isfinite(v).all() for v in (mem, prior, prior_grad, jacobian)):
        raise ValueError("registration coefficient arrays must be finite")
    if not isinstance(offset, (int, np.integer)) or offset < 0:
        raise ValueError("offset must be a non-negative integer")
    if not isinstance(subsample, (int, np.integer)) or subsample < 1:
        raise ValueError("subsample must be a positive integer")
    if offset >= min(mem.shape[:3]):
        raise ValueError("offset must be inside every spatial axis")
    if registered_classes is None:
        registered = np.ones(mem.shape[-1], dtype=bool)
    else:
        registered = np.asarray(registered_classes, dtype=bool)
        if registered.shape != (mem.shape[-1],):
            raise ValueError("registered_classes must have one value per class")

    cost = np.float32(0.0)
    norm = np.float32(0.0)
    gradient = np.zeros(6, dtype=np.float32)
    hessian = np.zeros(6, dtype=np.float32)
    norm_gradient = np.zeros(6, dtype=np.float32)
    norm_hessian = np.zeros(6, dtype=np.float32)
    for x in range(offset, mem.shape[0], subsample):
        for y in range(offset, mem.shape[1], subsample):
            for z in range(offset, mem.shape[2], subsample):
                for k in range(1, mem.shape[-1]):
                    p = prior[x, y, z, k]
                    if not registered[k] or p <= np.float32(0.0):
                        continue
                    m = mem[x, y, z, k]
                    m2 = np.float32(m * m)
                    p2 = np.float32(p * p)
                    cost = np.float32(cost + np.float32(np.float32(m2 * p) * p))
                    norm = np.float32(norm + p2)
                    dp = np.empty(6, dtype=np.float32)
                    for parameter in range(6):
                        value = np.float32(prior_grad[x, y, z, k, 0] * jacobian[x, y, z, 0, parameter])
                        value = np.float32(value + np.float32(prior_grad[x, y, z, k, 1] * jacobian[x, y, z, 1, parameter]))
                        dp[parameter] = np.float32(value + np.float32(prior_grad[x, y, z, k, 2] * jacobian[x, y, z, 2, parameter]))
                    for parameter in range(6):
                        gradient_term = np.float32(np.float32(2.0) * m)
                        gradient_term = np.float32(gradient_term * m)
                        gradient_term = np.float32(gradient_term * p)
                        gradient_term = np.float32(gradient_term * dp[parameter])
                        gradient[parameter] = np.float32(gradient[parameter] + gradient_term)
                        hessian_term = np.float32(np.float32(2.0) * m)
                        hessian_term = np.float32(hessian_term * m)
                        hessian_term = np.float32(hessian_term * dp[parameter])
                        hessian_term = np.float32(hessian_term * dp[parameter])
                        hessian[parameter] = np.float32(hessian[parameter] + hessian_term)
                        norm_gradient_term = np.float32(np.float32(2.0) * p)
                        norm_gradient_term = np.float32(norm_gradient_term * dp[parameter])
                        norm_gradient[parameter] = np.float32(norm_gradient[parameter] + norm_gradient_term)
                        norm_hessian_term = np.float32(np.float32(2.0) * dp[parameter])
                        norm_hessian_term = np.float32(norm_hessian_term * dp[parameter])
                        norm_hessian[parameter] = np.float32(norm_hessian[parameter] + norm_hessian_term)
    return RegistrationCoefficientSums(
        cost, norm, gradient, hessian, norm_gradient, norm_hessian
    )


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


def trilinear_prior_derivatives(
    shape_prior: ArrayLike,
    template_xyz: ArrayLike,
) -> NDArray[np.float32]:
    """Return source-compatible x/y/z derivatives of a trilinear atlas prior.

    This follows `ImageFunctions.linearInterpolation[XYZ]derivative`, which
    returns zero outside the closed domain `[0, size - 2]` on any axis. That
    boundary policy is deliberately different from the closest-edge value
    interpolation used for the registration energy.
    """
    values = np.asarray(shape_prior, dtype=np.float32)
    point = np.asarray(template_xyz, dtype=np.float32)
    if values.ndim != 3 or any(size < 2 for size in values.shape):
        raise ValueError("shape_prior must be a 3D array with dimensions >= 2")
    if point.shape != (3,) or not np.isfinite(point).all():
        raise ValueError("template_xyz must contain three finite coordinates")
    if not np.isfinite(values).all():
        raise ValueError("shape_prior must contain only finite values")
    if any(point[axis] < 0.0 or point[axis] > values.shape[axis] - 2 for axis in range(3)):
        return np.zeros(3, dtype=np.float32)

    xyz0 = tuple(int(np.floor(point[axis])) for axis in range(3))
    alpha, beta, gamma = (
        np.float32(point[axis] - np.float32(xyz0[axis])) for axis in range(3)
    )
    nalpha, nbeta, ngamma = (
        np.float32(np.float32(1.0) - t) for t in (alpha, beta, gamma)
    )
    x, y, z = xyz0

    def product(*factors: np.float32) -> np.float32:
        result = np.float32(factors[0])
        for factor in factors[1:]:
            result = np.float32(result * np.float32(factor))
        return result

    def difference(a: tuple[int, int, int], b: tuple[int, int, int]) -> np.float32:
        return np.float32(values[a] - values[b])

    dx_terms = (
        product(nbeta, ngamma, difference((x + 1, y, z), (x, y, z))),
        product(beta, ngamma, difference((x + 1, y + 1, z), (x, y + 1, z))),
        product(nbeta, gamma, difference((x + 1, y, z + 1), (x, y, z + 1))),
        product(beta, gamma, difference((x + 1, y + 1, z + 1), (x, y + 1, z + 1))),
    )
    dy_terms = (
        product(nalpha, ngamma, difference((x, y + 1, z), (x, y, z))),
        product(alpha, ngamma, difference((x + 1, y + 1, z), (x + 1, y, z))),
        product(nalpha, gamma, difference((x, y + 1, z + 1), (x, y, z + 1))),
        product(alpha, gamma, difference((x + 1, y + 1, z + 1), (x + 1, y, z + 1))),
    )
    dz_terms = (
        product(nalpha, nbeta, difference((x, y, z + 1), (x, y, z))),
        product(alpha, nbeta, difference((x + 1, y, z + 1), (x + 1, y, z))),
        product(nalpha, beta, difference((x, y + 1, z + 1), (x, y + 1, z))),
        product(alpha, beta, difference((x + 1, y + 1, z + 1), (x + 1, y + 1, z))),
    )
    result = np.zeros(3, dtype=np.float32)
    for axis, terms in enumerate((dx_terms, dy_terms, dz_terms)):
        total = np.float32(terms[0] + terms[1])
        total = np.float32(total + terms[2])
        result[axis] = np.float32(total + terms[3])
    return result


def rigid_image_to_template_jacobian(
    image_xyz: ArrayLike,
    transform_parameters: ArrayLike,
    image_center_xyz: ArrayLike,
    image_resolutions_mm: ArrayLike,
    template_resolutions_mm: ArrayLike,
    *,
    scale: float = 1.0,
) -> NDArray[np.float32]:
    """Return the source rigid-map coordinate Jacobian for one input voxel.

    Output axes are `(template coordinate, rigid parameter)` where parameters
    are inverse-Rodrigues x/y/z followed by x/y/z millimetre translations.
    The implementation follows the selected `RotationMatrix.derivatives` and
    `ParametricTransform.imageToTemplateDerivatives` float32 formulas.
    """
    voxel = np.asarray(image_xyz, dtype=np.float32)
    transform = np.asarray(transform_parameters, dtype=np.float32)
    center = np.asarray(image_center_xyz, dtype=np.float32)
    image_res = np.asarray(image_resolutions_mm, dtype=np.float32)
    template_res = np.asarray(template_resolutions_mm, dtype=np.float32)
    scale_value = np.float32(scale)
    if voxel.shape != (3,) or not np.isfinite(voxel).all():
        raise ValueError("image_xyz must contain three finite coordinates")
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("transform_parameters must contain six finite rigid parameters")
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("image_center_xyz must contain three finite values")
    if image_res.shape != (3,) or not np.isfinite(image_res).all() or np.any(image_res <= 0):
        raise ValueError("image_resolutions_mm must contain three positive finite values")
    if template_res.shape != (3,) or not np.isfinite(template_res).all() or np.any(template_res <= 0):
        raise ValueError("template_resolutions_mm must contain three positive finite values")
    if not np.isfinite(scale_value) or scale_value <= 0:
        raise ValueError("scale must be positive and finite")

    def f(value: float) -> np.float32:
        return np.float32(value)

    a, b, c = (f(v) for v in transform[:3])
    norm_square = f(f(f(a * a) + f(b * b)) + f(c * c))
    norm = f(np.sqrt(norm_square))
    if norm > f(1.0):
        magnitude = f(norm % f(2.0))
        if magnitude > f(1.0):
            a = f(f(f(magnitude * a) / norm) - f(f(f(2.0) * a) / norm))
            b = f(f(f(magnitude * b) / norm) - f(f(f(2.0) * b) / norm))
            c = f(f(f(magnitude * c) / norm) - f(f(f(2.0) * c) / norm))
        else:
            a, b, c = (f(f(v * magnitude) / norm) for v in (a, b, c))
    a2, b2, c2 = f(a * a), f(b * b), f(c * c)
    ab, ac, bc = f(a * b), f(a * c), f(b * c)
    x2 = f(f(f(1.0) - a2) - b2)
    x2 = f(x2 - c2)
    x = f(np.sqrt(x2))
    if x == f(0.0):
        raise ValueError("rigid rotation derivative is singular at unit parameter norm")
    ax, bx, cx = f(x * a), f(x * b), f(x * c)

    def rotation_derivative(da: np.float32, db: np.float32, dc: np.float32) -> NDArray[np.float32]:
        out = np.zeros((3, 3), dtype=np.float32)
        out[0, 0] = f(f(-4.0) * f(f(b * db) + f(c * dc)))
        q = f(f(f(bx - ac) * da) + f(f(ax - bc) * db))
        q = f(q + f(f(x2 - c2) * dc))
        out[1, 0] = f(f(f(2.0) / x) * q)
        q = f(f(f(cx + ab) * da) + f(f(b2 - x2) * db))
        q = f(q + f(f(ax + bc) * dc))
        out[2, 0] = f(f(f(2.0) / x) * q)
        q = f(f(f(bx + ac) * da) + f(f(ax + bc) * db))
        q = f(q + f(f(c2 - x2) * dc))
        out[0, 1] = f(f(f(2.0) / x) * q)
        out[1, 1] = f(f(f(-4.0) * f(f(a * da) + f(c * dc))))
        q = f(f(f(x2 - a2) * da) + f(f(cx - ab) * db))
        q = f(q + f(f(bx - ac) * dc))
        out[2, 1] = f(f(f(2.0) / x) * q)
        q = f(f(f(cx - ab) * da) + f(f(x2 - b2) * db))
        q = f(q + f(f(ax - bc) * dc))
        out[0, 2] = f(f(f(2.0) / x) * q)
        q = f(f(f(a2 - x2) * da) + f(f(cx + ab) * db))
        q = f(q + f(f(bx + ac) * dc))
        out[1, 2] = f(f(f(2.0) / x) * q)
        out[2, 2] = f(f(f(-4.0) * f(f(a * da) + f(b * db))))
        return out

    derivatives = (
        rotation_derivative(f(1.0), f(0.0), f(0.0)),
        rotation_derivative(f(0.0), f(1.0), f(0.0)),
        rotation_derivative(f(0.0), f(0.0), f(1.0)),
    )
    physical = tuple(
        f(f(f(voxel[axis] * scale_value) - center[axis]) * image_res[axis])
        for axis in range(3)
    )
    jacobian = np.zeros((3, 6), dtype=np.float32)
    for parameter, derivative in enumerate(derivatives):
        for axis in range(3):
            value = f(derivative[axis, 0] * physical[0])
            value = f(value + f(derivative[axis, 1] * physical[1]))
            value = f(value + f(derivative[axis, 2] * physical[2]))
            value = f(value / template_res[axis])
            jacobian[axis, parameter] = f(value / scale_value)
    jacobian[0, 3] = f(f(f(1.0) / template_res[0]) / scale_value)
    jacobian[1, 4] = f(f(f(1.0) / template_res[1]) / scale_value)
    jacobian[2, 5] = f(f(f(1.0) / template_res[2]) / scale_value)
    return jacobian


def synthetic_registration_objective_inputs() -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Return the deterministic 5×4×3×4 Java-reference objective phantom."""
    x, y, z, k = np.indices((5, 4, 3, 4), dtype=np.int32)
    memberships = (((7 * x + 3 * y + 5 * z + 2 * k) % 19) / np.float32(18.0)).astype(np.float32)
    priors = (((2 * x + 5 * y + 3 * z + 4 * k) % 11) / np.float32(10.0)).astype(np.float32)
    return memberships, priors
