"""Selected MaCRUISE ACE stages ported from the supplied VandyCRUISE code.

Source anchors include Blake Lucas's VandyCRUISE `AnisotropicDiffusion`,
`Morphology`, and `CleanACE`, and Yuankai Huo's MaACE workflow adapter.
Upstream notices are preserved under `licenses/upstream/`.
"""

from __future__ import annotations

from dataclasses import dataclass
import heapq

import nibabel as nib
import numpy as np
from scipy import ndimage

from pymacruise.image_io import normalize_image_to_ras

_STRUCTURE_18 = ndimage.generate_binary_structure(3, 2)
_STRUCTURE_6 = ndimage.generate_binary_structure(3, 1)


@dataclass(frozen=True)
class ACESeedFields:
    """Binary and CSF initialization products before the ACE PDE stages."""

    initial_csf_speed: np.ndarray
    prohibited_region: np.ndarray
    closed_combined_wm: np.ndarray
    naive_csf: np.ndarray


@dataclass(frozen=True)
class ACEArrayResult:
    """ACE outputs on one already-normalized 3D voxel grid.

    ``skeleton`` is the cleaned ACE CSF membership returned before the
    separate MaCRUISE thinning step (roadmap P27).
    """

    enhanced_gray_matter: np.ndarray
    skeleton: np.ndarray


@dataclass(frozen=True)
class ACEImageResult:
    """RAS-oriented NIfTI ACE outputs on the input physical grid."""

    enhanced_gray_matter: nib.Nifti1Image
    skeleton: nib.Nifti1Image


def ace_seed_fields(
    gray_matter: np.ndarray,
    white_matter: np.ndarray,
    signed_distance: np.ndarray,
    *,
    intensity_threshold: float = 126.9,
) -> ACESeedFields:
    """Port MaACE's membership thresholds and initial CSF speed mask."""
    gray = np.asarray(gray_matter, dtype=np.float32)
    white = np.asarray(white_matter, dtype=np.float32)
    phi = np.asarray(signed_distance, dtype=np.float32)
    if gray.ndim != 3 or gray.shape != white.shape or gray.shape != phi.shape:
        raise ValueError("ACE gray, white and signed-distance arrays must share a 3D grid")
    if not np.isfinite(gray).all() or not np.isfinite(white).all() or not np.isfinite(phi).all():
        raise ValueError("ACE seed inputs must contain only finite values")
    if not np.isfinite(intensity_threshold):
        raise ValueError("intensity_threshold must be finite")

    initial_mask = white >= 30.0
    combined = white + gray
    combined_wm = (combined >= 210.0).astype(np.uint8)
    truncated_sum = np.trunc(combined).astype(np.int32)
    naive_csf = np.where(combined > 254.0, 0, 254 - truncated_sum).astype(np.int32)
    closed_wm = ace_close_wm(combined_wm, iterations=4)
    candidate = (closed_wm == 1) & (combined_wm == 0) & (phi > 0.0)
    initial_csf = np.zeros(gray.shape, dtype=np.float32)
    initial_csf[candidate] = np.sqrt(naive_csf[candidate] / 254.0).astype(np.float32)
    prohibited = (phi <= 0.0) | initial_mask
    return ACESeedFields(initial_csf, prohibited, closed_wm, naive_csf)


def ace_anisotropic_diffusion(
    initial: np.ndarray,
    guide: np.ndarray,
    *,
    tensor_alpha: float = 1.6,
    iterations: int = 25,
) -> np.ndarray:
    """Port ACE's structure-tensor diffusion and monotone update rule."""
    values = np.asarray(initial, dtype=np.float32)
    guide_values = np.asarray(guide, dtype=np.float32)
    if values.ndim != 3 or values.shape != guide_values.shape:
        raise ValueError("ACE diffusion inputs must share one 3D grid")
    if not np.isfinite(values).all() or not np.isfinite(guide_values).all():
        raise ValueError("ACE diffusion inputs must be finite")
    if not np.isfinite(tensor_alpha) or tensor_alpha <= 0:
        raise ValueError("tensor_alpha must be finite and positive")
    if iterations < 0:
        raise ValueError("iterations must be non-negative")

    gradient = _central_gradient_mirror(guide_values)
    px, py, pz = gradient
    tensor = [
        ace_gaussian_blur(px * px, tensor_alpha),
        ace_gaussian_blur(px * py, tensor_alpha),
        ace_gaussian_blur(px * pz, tensor_alpha),
        ace_gaussian_blur(py * py, tensor_alpha),
        ace_gaussian_blur(py * pz, tensor_alpha),
        ace_gaussian_blur(pz * pz, tensor_alpha),
    ]
    a11, a12, a13, a22, a23, a33 = tensor
    vector = _principal_tensor_vector(a11, a12, a13, a22, a23, a33)
    vx, vy, vz = vector
    diffusion = (
        np.float32(0.999) - np.float32(0.998) * vx * vx,
        -np.float32(0.998) * vx * vy,
        -np.float32(0.998) * vx * vz,
        np.float32(0.999) - np.float32(0.998) * vy * vy,
        -np.float32(0.998) * vy * vz,
        np.float32(0.999) - np.float32(0.998) * vz * vz,
    )
    d11, d12, d13, d22, d23, d33 = diffusion

    output = values.copy()
    for _ in range(iterations):
        backward = _backward_gradient_mirror(output)
        gx, gy, gz = backward
        flux_x = d11 * gx + d12 * gy + d13 * gz
        flux_y = d12 * gx + d22 * gy + d23 * gz
        flux_z = d13 * gx + d23 * gy + d33 * gz
        divergence = (
            _forward_difference_mirror(flux_x, 0)
            + _forward_difference_mirror(flux_y, 1)
            + _forward_difference_mirror(flux_z, 2)
        )
        updated = output + np.float32(0.1) * divergence
        output = np.clip(np.maximum(updated, output), 0.0, 1.0).astype(np.float32)
    return output


def ace_update_speed(
    diffused_csf: np.ndarray,
    *,
    multi_atlas_csf: np.ndarray | None = None,
    gaussian_alpha: float = 1.8,
) -> np.ndarray:
    """Apply V3 after-diffusion CSF merge and ACE speed-function blur."""
    csf = np.asarray(diffused_csf, dtype=np.float32).copy()
    if csf.ndim != 3 or not np.isfinite(csf).all():
        raise ValueError("diffused_csf must be a finite 3D array")
    if multi_atlas_csf is not None:
        updated = np.asarray(multi_atlas_csf)
        if updated.shape != csf.shape:
            raise ValueError("multi-atlas CSF must share the diffused CSF grid")
        if not np.isfinite(updated).all():
            raise ValueError("multi-atlas CSF must be finite")
        # The selected wrapper uses Merge (option 0) after diffusion: any
        # positive atlas CSF sets the initial speed function to one.
        csf[updated > 0] = 1.0
    speed = 1.0 - 0.9 * csf
    return ace_gaussian_blur(speed, gaussian_alpha)


def ace_gaussian_blur(image: np.ndarray, alpha: float) -> np.ndarray:
    """Run the ACE Deriche recursive Gaussian with source ``alpha`` semantics.

    The source documents ``alpha = 1.4 / sigma`` and applies its recursive
    filter along voxel axes 2, 1, then 0. This is intentionally distinct from
    SciPy's sampled-kernel `gaussian_filter`.
    """
    values = np.asarray(image, dtype=np.float32)
    if values.ndim != 3:
        raise ValueError("ACE Gaussian blur requires a 3D scalar array")
    if any(length < 2 for length in values.shape):
        raise ValueError("ACE Gaussian blur requires each dimension to be >= 2")
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be a finite positive value")

    d1 = float(np.exp(-alpha))
    d2 = d1 * d1
    b0 = (1.0 - d1) ** 2 / (1.0 + 2.0 * alpha * d1 - d2)
    b1 = 2.0 * d1
    b2 = -d2
    a1 = b0
    a2 = b0 * d1 * (alpha - 1.0)
    a3 = b0 * d1 * (alpha + 1.0)
    a4 = -b0 * d2

    output = values
    for axis in (2, 1, 0):
        moved = np.moveaxis(output, axis, -1)
        forward = np.empty_like(moved)
        backward = np.zeros_like(moved)
        forward[..., 0] = (a1 * moved[..., 0]).astype(np.float32)
        forward[..., 1] = (
            a1 * moved[..., 1] + a2 * moved[..., 0] + b1 * forward[..., 0]
        ).astype(np.float32)
        for index in range(2, moved.shape[-1]):
            forward[..., index] = (
                a1 * moved[..., index]
                + a2 * moved[..., index - 1]
                + b1 * forward[..., index - 1]
                + b2 * forward[..., index - 2]
            ).astype(np.float32)

        backward[..., -1] = 0
        backward[..., -2] = (
            a3 * moved[..., -1] + b1 * backward[..., -1]
        ).astype(np.float32)
        for index in range(moved.shape[-1] - 3, -1, -1):
            backward[..., index] = (
                a3 * moved[..., index + 1]
                + a4 * moved[..., index + 2]
                + b1 * backward[..., index + 1]
                + b2 * backward[..., index + 2]
            ).astype(np.float32)
        output = np.moveaxis(forward + backward, -1, axis)
    return np.asarray(output, dtype=np.float32)


def ace_close_wm(wm_membership: np.ndarray, iterations: int = 4) -> np.ndarray:
    """Reproduce ACE's binary closing with the JIST 18-neighbor footprint."""
    values = np.asarray(wm_membership)
    if values.ndim != 3:
        raise ValueError("ACE WM morphology requires a 3D array")
    if iterations < 0:
        raise ValueError("iterations must be non-negative")
    # The JIST adapter converts to ImageDataInt before applying its binary
    # morphology, so foreground membership is integer-nonzero.
    foreground = values.astype(np.int32) != 0
    for _ in range(iterations):
        foreground = ndimage.binary_dilation(
            foreground, structure=_STRUCTURE_18, border_value=0
        )
    for _ in range(iterations):
        # The Java implementation ignores out-of-image neighbors during
        # erosion, equivalent to a foreground border value of one.
        foreground = ndimage.binary_erosion(
            foreground, structure=_STRUCTURE_18, border_value=1
        )
    return foreground.astype(np.uint8)


def clean_ace_components(
    csf_membership: np.ndarray,
    *,
    minimum_size_exclusive: int = 60,
) -> np.ndarray:
    """Remove ACE components with size <= 60 using source 18-connectivity.

    The source scans seed voxels only in the interior, while a seeded component
    may extend to the image border. Border-only components are left untouched.
    """
    values = np.asarray(csf_membership)
    if values.ndim != 3:
        raise ValueError("ACE component cleanup requires a 3D array")
    if minimum_size_exclusive < 0:
        raise ValueError("minimum_size_exclusive must be non-negative")
    foreground = values > 0
    labels, count = ndimage.label(foreground, structure=_STRUCTURE_18)
    if count == 0:
        return values.copy()
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    interior_labels = np.unique(labels[1:-1, 1:-1, 1:-1])
    remove_ids = interior_labels[
        (interior_labels > 0) & (sizes[interior_labels] <= minimum_size_exclusive)
    ]
    if remove_ids.size == 0:
        return values.copy()
    result = values.copy()
    result[np.isin(labels, remove_ids)] = 0
    return result


def ace_enhance_arrays(
    gray_matter: np.ndarray,
    white_matter: np.ndarray,
    *,
    multi_atlas_csf: np.ndarray | None = None,
    intensity_threshold: float = 126.9,
    tensor_alpha: float = 1.6,
    gaussian_alpha: float = 1.8,
    diffusion_iterations: int = 25,
) -> ACEArrayResult:
    """Run the selected V3 MaACE enhancement path on aligned 3D arrays.

    Inputs must already use the same RAS voxel grid; the NIfTI-facing wrapper
    performs orientation normalization and emits output images in a later API
    layer. Updated multi-atlas CSF follows the selected V3 after-diffusion
    merge setting.
    """
    gray = np.asarray(gray_matter, dtype=np.float32)
    white = np.asarray(white_matter, dtype=np.float32)
    if gray.ndim != 3 or gray.shape != white.shape:
        raise ValueError("ACE gray and white matter arrays must share one 3D grid")
    if not np.isfinite(gray).all() or not np.isfinite(white).all():
        raise ValueError("ACE tissue inputs must be finite")
    if any(length < 2 for length in gray.shape):
        raise ValueError("ACE inputs require dimensions >= 2")
    if not np.isfinite(intensity_threshold):
        raise ValueError("intensity_threshold must be finite")
    if not np.isfinite(tensor_alpha) or tensor_alpha <= 0:
        raise ValueError("tensor_alpha must be finite and positive")
    if not np.isfinite(gaussian_alpha) or gaussian_alpha <= 0:
        raise ValueError("gaussian_alpha must be finite and positive")
    if diffusion_iterations < 0:
        raise ValueError("diffusion_iterations must be non-negative")
    if multi_atlas_csf is not None:
        atlas_csf = np.asarray(multi_atlas_csf)
        if atlas_csf.shape != gray.shape or not np.isfinite(atlas_csf).all():
            raise ValueError("multi_atlas_csf must be finite and share the tissue grid")
    else:
        atlas_csf = None

    initial_phi = _signed_distance_fmm(intensity_threshold - white, 1000.0)
    seeds = ace_seed_fields(gray, white, initial_phi, intensity_threshold=intensity_threshold)
    guide = ace_gaussian_blur(initial_phi, gaussian_alpha)
    diffused_csf = ace_anisotropic_diffusion(
        seeds.initial_csf_speed,
        guide,
        tensor_alpha=tensor_alpha,
        iterations=diffusion_iterations,
    )
    speed = ace_update_speed(
        diffused_csf,
        multi_atlas_csf=atlas_csf,
        gaussian_alpha=gaussian_alpha,
    )
    evolved_phi = narrow_band_ace(
        initial_phi,
        delta_t=0.1,
        max_iterations=1,
        inner_steps=10,
        curvature_threshold=5.0,
    )
    weighted_distance, _ = _distance_and_characteristic_fmm(evolved_phi, speed)

    dx, dy, dz = _central_gradient_mirror(weighted_distance)
    magnitude = np.sqrt(dx * dx + dy * dy + dz * dz) * speed
    eligible = (weighted_distance > 0.25) & ~seeds.prohibited_region
    memberships = np.zeros(gray.shape, dtype=np.int32)
    selected = eligible & (magnitude < 0.85)
    memberships[selected] = np.trunc(254.0 * (magnitude[selected] + 0.05)).astype(np.int32)
    memberships = clean_ace_components(memberships, minimum_size_exclusive=60)
    memberships[memberships > 0] = 254

    enhanced = gray.copy()
    enhanced[memberships > 0] = 0.0
    border = np.zeros(gray.shape, dtype=bool)
    border[[0, -1], :, :] = True
    border[:, [0, -1], :] = True
    border[:, :, [0, -1]] = True
    enhanced[border] = 0.0
    memberships[border] = 0
    return ACEArrayResult(enhanced, memberships)


def ace_enhance(
    gray_matter: nib.spatialimages.SpatialImage,
    white_matter: nib.spatialimages.SpatialImage,
    *,
    multi_atlas_csf: nib.spatialimages.SpatialImage | None = None,
    intensity_threshold: float = 126.9,
    tensor_alpha: float = 1.6,
    gaussian_alpha: float = 1.8,
    diffusion_iterations: int = 25,
) -> ACEImageResult:
    """Run ACE on NIfTI inputs, reorienting to RAS and warning when needed."""
    gray_image = _as_ras_3d(gray_matter, "ACE gray matter")
    white_image = _as_ras_3d(white_matter, "ACE white matter")
    _same_ras_grid(gray_image, white_image, "ACE gray matter", "ACE white matter")
    atlas_image = None
    if multi_atlas_csf is not None:
        atlas_image = _as_ras_3d(multi_atlas_csf, "ACE multi-atlas CSF")
        _same_ras_grid(gray_image, atlas_image, "ACE gray matter", "ACE multi-atlas CSF")

    arrays = ace_enhance_arrays(
        np.asarray(gray_image.dataobj),
        np.asarray(white_image.dataobj),
        multi_atlas_csf=None if atlas_image is None else np.asarray(atlas_image.dataobj),
        intensity_threshold=intensity_threshold,
        tensor_alpha=tensor_alpha,
        gaussian_alpha=gaussian_alpha,
        diffusion_iterations=diffusion_iterations,
    )
    return ACEImageResult(
        _make_output(arrays.enhanced_gray_matter, gray_image),
        _make_output(arrays.skeleton, gray_image),
    )


def _central_gradient_mirror(values: np.ndarray) -> tuple[np.ndarray, ...]:
    """Source centered difference with edge samples replicated at the boundary."""
    gradient = []
    for axis in range(3):
        before = np.maximum(np.arange(values.shape[axis]) - 1, 0)
        after = np.minimum(np.arange(values.shape[axis]) + 1, values.shape[axis] - 1)
        gradient.append(
            (np.take(values, after, axis=axis) - np.take(values, before, axis=axis))
            * np.float32(0.5)
        )
    return tuple(np.asarray(item, dtype=np.float32) for item in gradient)


def _backward_gradient_mirror(values: np.ndarray) -> tuple[np.ndarray, ...]:
    """Source backward difference; the first voxel has a zero difference."""
    gradient = []
    for axis in range(3):
        previous = np.maximum(np.arange(values.shape[axis]) - 1, 0)
        gradient.append(values - np.take(values, previous, axis=axis))
    return tuple(np.asarray(item, dtype=np.float32) for item in gradient)


def _forward_difference_mirror(values: np.ndarray, axis: int) -> np.ndarray:
    """Source forward difference; the last voxel has a zero difference."""
    following = np.minimum(np.arange(values.shape[axis]) + 1, values.shape[axis] - 1)
    return np.asarray(np.take(values, following, axis=axis) - values, dtype=np.float32)


def _principal_tensor_vector(
    a11: np.ndarray,
    a12: np.ndarray,
    a13: np.ndarray,
    a22: np.ndarray,
    a23: np.ndarray,
    a33: np.ndarray,
) -> tuple[np.ndarray, ...]:
    """Select the largest-eigenvalue vector using the source's 10 power steps."""
    # Match the Java branch and tie behavior exactly: ties involving a33 favor
    # z, while an a11/a22 tie favors whichever strictly exceeds a33.
    px = np.zeros_like(a11, dtype=np.float32)
    py = np.zeros_like(a11, dtype=np.float32)
    pz = np.zeros_like(a11, dtype=np.float32)
    choose_x = (a11 > a22) & (a11 > a33)
    choose_y = (a11 <= a22) & (a22 > a33)
    choose_z = ~(choose_x | choose_y)
    px[choose_x] = 1.0
    py[choose_y] = 1.0
    pz[choose_z] = 1.0

    for _ in range(10):
        nx = a11 * px + a12 * py + a13 * pz
        ny = a12 * px + a22 * py + a23 * pz
        nz = a13 * px + a23 * py + a33 * pz
        magnitude = np.sqrt(nx * nx + ny * ny + nz * nz).astype(np.float32)
        nonzero = magnitude > 0
        # The Java implementation divides by zero for a zero tensor. Preserve
        # its initialized basis vector there, avoiding NaNs in otherwise valid
        # uniform-guide volumes.
        px = np.divide(nx, magnitude, out=px.copy(), where=nonzero)
        py = np.divide(ny, magnitude, out=py.copy(), where=nonzero)
        pz = np.divide(nz, magnitude, out=pz.copy(), where=nonzero)
    return px, py, pz


def _signed_distance_fmm(level_set: np.ndarray, max_distance: float = 1000.0) -> np.ndarray:
    """Port VandyCRUISE's six-neighbor, unit-speed signed fast-marching field.

    Zero voxels and sign-crossing voxels form the initial accepted set. Crossing
    distances use linear subvoxel interpolation; the Eikonal update then marches
    independently outward from both sides of the interface.
    """
    phi = np.asarray(level_set, dtype=np.float32)
    if phi.ndim != 3 or any(size < 2 for size in phi.shape):
        raise ValueError("ACE signed distance requires a 3D array with dimensions >= 2")
    if not np.isfinite(phi).all():
        raise ValueError("ACE signed-distance input must be finite")
    if not np.isfinite(max_distance) or max_distance <= 0:
        raise ValueError("max_distance must be finite and positive")

    alive = phi == 0
    distance = np.full(phi.shape, np.inf, dtype=np.float64)
    distance[alive] = 0.0
    seeds = np.argwhere(~alive)
    for point in seeds:
        index = tuple(point)
        value = float(phi[index])
        offsets = _source_crossing_offsets(phi, index)
        if offsets:
            distance[index] = 1.0 / np.sqrt(sum(1.0 / (offset * offset) for offset in offsets))
            alive[index] = True

    state = np.zeros(phi.shape, dtype=np.uint8)  # 0 far, 1 narrow band, 2 accepted
    state[alive] = 2
    heap: list[tuple[float, int, int, int]] = []

    def update(candidate: tuple[int, int, int]) -> float:
        terms = []
        for axis in range(3):
            neighbor_values = []
            for delta in (-1, 1):
                neighbor_index = list(candidate)
                neighbor_index[axis] += delta
                if 0 <= neighbor_index[axis] < phi.shape[axis]:
                    ni = tuple(neighbor_index)
                    if state[ni] == 2:
                        neighbor_values.append(distance[ni])
            if neighbor_values:
                terms.append(min(neighbor_values))
        if not terms:
            return float("inf")
        terms.sort()
        total = 0.0
        square_total = 0.0
        root = terms[0] + 1.0
        for count, term in enumerate(terms, start=1):
            total += term
            square_total += term * term
            discriminant = total * total - count * (square_total - 1.0)
            root = (total + np.sqrt(max(0.0, discriminant))) / count
            if count == len(terms) or root <= terms[count]:
                break
        return float(root)

    for point in np.argwhere(alive):
        current = tuple(int(v) for v in point)
        for axis in range(3):
            for delta in (-1, 1):
                neighbor = list(current)
                neighbor[axis] += delta
                if not 0 <= neighbor[axis] < phi.shape[axis]:
                    continue
                candidate = tuple(neighbor)
                if state[candidate] == 2:
                    continue
                tentative = update(candidate)
                if tentative < distance[candidate]:
                    distance[candidate] = tentative
                    state[candidate] = 1
                    heapq.heappush(heap, (tentative, *candidate))

    while heap:
        value, i, j, k = heapq.heappop(heap)
        current = (i, j, k)
        if state[current] == 2 or value != distance[current]:
            continue
        if value > max_distance:
            break
        state[current] = 2
        for axis in range(3):
            for delta in (-1, 1):
                neighbor = [i, j, k]
                neighbor[axis] += delta
                if not 0 <= neighbor[axis] < phi.shape[axis]:
                    continue
                candidate = tuple(neighbor)
                if state[candidate] == 2:
                    continue
                tentative = update(candidate)
                if tentative < distance[candidate]:
                    distance[candidate] = tentative
                    state[candidate] = 1
                    heapq.heappush(heap, (tentative, *candidate))

    distance[state != 2] = max_distance
    signed = np.where(phi < 0, -distance, distance)
    return np.asarray(signed, dtype=np.float32)


def _shift_clamped(values: np.ndarray, axis: int, offset: int) -> np.ndarray:
    indices = np.clip(
        np.arange(values.shape[axis]) + offset, 0, values.shape[axis] - 1
    )
    return np.take(values, indices, axis=axis)


def _source_crossing_offsets(values: np.ndarray, point: tuple[int, int, int]) -> list[float]:
    """Return one source-compatible subvoxel crossing per coordinate axis."""
    center = float(values[point])
    offsets = []
    for axis in range(3):
        opposite = []
        for delta in (-1, 1):
            neighbor = list(point)
            neighbor[axis] += delta
            if 0 <= neighbor[axis] < values.shape[axis]:
                adjacent = float(values[tuple(neighbor)])
                if center * adjacent < 0:
                    opposite.append(adjacent)
        if opposite:
            # The Java code resolves a two-sided crossing to the neighbor with
            # the greater absolute level-set value (ties choose +axis).
            adjacent = opposite[0] if len(opposite) == 1 or abs(opposite[0]) > abs(opposite[1]) else opposite[1]
            offsets.append(abs(center / (center - adjacent)))
    return offsets


def narrow_band_ace(
    signed_distance: np.ndarray,
    *,
    delta_t: float = 0.1,
    max_iterations: int = 1,
    inner_steps: int = 10,
    curvature_threshold: float = 5.0,
) -> np.ndarray:
    """Port VandyCRUISE NarrowBandACE and its fast-marching reinitialization."""
    phi = np.asarray(signed_distance, dtype=np.float32).copy()
    if phi.ndim != 3 or any(size < 2 for size in phi.shape):
        raise ValueError("NarrowBandACE requires a 3D array with dimensions >= 2")
    if not np.isfinite(phi).all():
        raise ValueError("NarrowBandACE input must be finite")
    if not np.isfinite(delta_t) or delta_t < 0:
        raise ValueError("delta_t must be finite and non-negative")
    if max_iterations < 0 or inner_steps < 0:
        raise ValueError("iteration counts must be non-negative")
    if not np.isfinite(curvature_threshold):
        raise ValueError("curvature_threshold must be finite")

    threshold = abs(float(curvature_threshold))
    for _ in range(max_iterations):
        band = np.abs(phi) <= 6.0
        mines = band & (np.abs(phi) > 4.0)
        mine_signs = np.sign(phi[mines])
        coordinates = np.nonzero(band)
        for _ in range(inner_steps):
            # Source mirror-boundary first/second/central differences. Axes
            # correspond to Java's (row, column, slice) as (z, y, x).
            west, east = _shift_clamped(phi, 2, -1), _shift_clamped(phi, 2, 1)
            north, south = _shift_clamped(phi, 1, -1), _shift_clamped(phi, 1, 1)
            back, front = _shift_clamped(phi, 0, -1), _shift_clamped(phi, 0, 1)
            dmx, dpx = phi - west, east - phi
            dmy, dpy = phi - north, south - phi
            dmz, dpz = phi - back, front - phi
            d0x, d0y, d0z = (east - west) / 2, (south - north) / 2, (front - back) / 2
            dxx, dyy, dzz = west + east - 2 * phi, north + south - 2 * phi, back + front - 2 * phi
            dxy = (
                _shift_clamped(_shift_clamped(phi, 1, 1), 2, 1)
                + _shift_clamped(_shift_clamped(phi, 1, -1), 2, -1)
                - _shift_clamped(_shift_clamped(phi, 1, 1), 2, -1)
                - _shift_clamped(_shift_clamped(phi, 1, -1), 2, 1)
            ) / 4
            dxz = (
                _shift_clamped(_shift_clamped(phi, 0, 1), 2, 1)
                + _shift_clamped(_shift_clamped(phi, 0, -1), 2, -1)
                - _shift_clamped(_shift_clamped(phi, 0, 1), 2, -1)
                - _shift_clamped(_shift_clamped(phi, 0, -1), 2, 1)
            ) / 4
            dyz = (
                _shift_clamped(_shift_clamped(phi, 0, 1), 1, 1)
                + _shift_clamped(_shift_clamped(phi, 0, -1), 1, -1)
                - _shift_clamped(_shift_clamped(phi, 0, 1), 1, -1)
                - _shift_clamped(_shift_clamped(phi, 0, -1), 1, 1)
            ) / 4
            grad_sq = d0x * d0x + d0y * d0y + d0z * d0z
            grad_norm = np.sqrt(grad_sq) + 0.0001
            curvature = (
                (dyy + dzz) * d0x * d0x
                + (dxx + dzz) * d0y * d0y
                + (dxx + dyy) * d0z * d0z
                - 2 * (d0x * d0y * dxy + d0x * d0z * dxz + d0y * d0z * dyz)
            ) / (grad_norm * grad_norm * grad_norm)
            curvature = np.clip(curvature, -threshold, threshold)
            delta_positive = (
                np.maximum(dmx, 0) ** 2
                + np.minimum(dpx, 0) ** 2
                + np.maximum(dmy, 0) ** 2
                + np.minimum(dpy, 0) ** 2
                + np.maximum(dmz, 0) ** 2
                + np.minimum(dpz, 0) ** 2
            )
            value = np.sqrt(delta_positive) * 0.01 - 0.95 * curvature * grad_norm
            updated = phi.copy()
            updated[coordinates] = phi[coordinates] - delta_t * value[coordinates]
            phi = updated
            if mines[tuple(coordinates)].size and np.any(np.sign(phi[mines]) != mine_signs):
                break
        phi = _signed_distance_fmm(phi, max_distance=4.0)
    return phi


def _distance_and_characteristic_fmm(
    level_set: np.ndarray, speed: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Port VandyCRUISE's speed-weighted distance/characteristic fast march."""
    phi = np.asarray(level_set, dtype=np.float32)
    velocity = np.asarray(speed, dtype=np.float32)
    if phi.ndim != 3 or phi.shape != velocity.shape:
        raise ValueError("weighted fast-marching inputs must share one 3D grid")
    if not np.isfinite(phi).all() or not np.isfinite(velocity).all() or np.any(velocity < 0):
        raise ValueError("weighted fast-marching inputs must be finite with non-negative speed")
    if any(size < 2 for size in phi.shape):
        raise ValueError("weighted fast marching requires dimensions >= 2")

    shape = phi.shape
    time = np.full(shape, np.inf, dtype=np.float64)
    characteristic = np.zeros(shape, dtype=np.float64)
    state = np.zeros(shape, dtype=np.uint8)  # far, narrow band, accepted
    for point in np.ndindex(shape):
        value = float(phi[point])
        if value == 0:
            time[point] = 0.0
            state[point] = 2
            continue
        offsets = _source_crossing_offsets(phi, point)
        if offsets:
            char = 1.0 / np.sqrt(sum(1.0 / (offset * offset) for offset in offsets))
            characteristic[point] = char
            time[point] = char / velocity[point] if velocity[point] > 0 else 10000.0
            state[point] = 2

    if not np.any(state == 2):
        # The Java routine leaves the default all-zero distance and
        # characteristic arrays unchanged when the input has no interface.
        zeros = np.zeros(shape, dtype=np.float32)
        return zeros, zeros.copy()

    def march(point: tuple[int, int, int]) -> tuple[float, float]:
        directions: list[tuple[tuple[float, float], ...]] = []
        min_time = float("inf")
        min_char = 0.0
        for axis in range(3):
            pair = []
            for delta in (-1, 1):
                neighbor = list(point)
                neighbor[axis] += delta
                if 0 <= neighbor[axis] < shape[axis]:
                    ni = tuple(neighbor)
                    if state[ni] == 2:
                        pair.append((time[ni], characteristic[ni]))
                else:
                    pair.append((float("inf"), 0.0))
            finite_alive = [item for item in pair if np.isfinite(item[0]) or item[0] == 10000.0]
            if finite_alive:
                selected = min(finite_alive, key=lambda item: item[0])
                if selected[0] < min_time:
                    min_time, min_char = selected
            directions.append(tuple(pair))

        local_speed = float(velocity[point])
        if local_speed > 0:
            min_time += 1.0 / local_speed
        summation = 0.0
        squares = 0.0
        char_sum = 0.0
        char_time_sum = 0.0
        count = 0
        for pair in directions:
            finite_alive = [item for item in pair if np.isfinite(item[0]) or item[0] == 10000.0]
            if not finite_alive:
                continue
            chosen = min(finite_alive, key=lambda item: item[0])
            if chosen[0] <= min_time:
                summation += chosen[0]
                squares += chosen[0] * chosen[0]
                char_sum += chosen[1]
                char_time_sum += chosen[1] * chosen[0]
                count += 1
        if count == 0:
            return min_time, min_char + 1.0
        if local_speed == 0:
            return 10000.0, char_sum / count + 1.0 / np.sqrt(count)
        inverse_speed = 1.0 / local_speed
        discriminant = summation * summation - count * (squares - inverse_speed * inverse_speed)
        if discriminant < 0:
            return min_time, min_char + 1.0
        result = (summation + np.sqrt(discriminant)) / count
        denominator = count * result - summation
        char = (inverse_speed + char_sum * result - char_time_sum) / denominator
        return float(result), float(char)

    heap: list[tuple[float, tuple[int, int, int]]] = []
    for point in np.argwhere(state == 2):
        accepted = tuple(int(value) for value in point)
        for axis in range(3):
            for delta in (-1, 1):
                neighbor = list(accepted)
                neighbor[axis] += delta
                if not 0 <= neighbor[axis] < shape[axis]:
                    continue
                candidate = tuple(neighbor)
                if state[candidate] == 2:
                    continue
                value, char = march(candidate)
                if value < time[candidate]:
                    time[candidate], characteristic[candidate] = value, char
                    state[candidate] = 1
                    heapq.heappush(heap, (value, candidate))

    while heap:
        value, point = heapq.heappop(heap)
        if state[point] == 2 or value != time[point]:
            continue
        state[point] = 2
        for axis in range(3):
            for delta in (-1, 1):
                neighbor = list(point)
                neighbor[axis] += delta
                if not 0 <= neighbor[axis] < shape[axis]:
                    continue
                candidate = tuple(neighbor)
                if state[candidate] == 2:
                    continue
                tentative, char = march(candidate)
                if tentative < time[candidate]:
                    time[candidate], characteristic[candidate] = tentative, char
                    state[candidate] = 1
                    heapq.heappush(heap, (tentative, candidate))

    time[phi < 0] *= -1
    return np.asarray(time, dtype=np.float32), np.asarray(characteristic, dtype=np.float32)


def _same_ras_grid(
    reference: nib.spatialimages.SpatialImage,
    candidate: nib.spatialimages.SpatialImage,
    reference_name: str,
    candidate_name: str,
) -> None:
    if reference.shape != candidate.shape or not np.allclose(
        reference.affine, candidate.affine, rtol=0.0, atol=1e-5
    ):
        raise ValueError(f"{reference_name} and {candidate_name} must share one physical grid")


def _as_ras_3d(
    image: nib.spatialimages.SpatialImage,
    name: str,
) -> nib.spatialimages.SpatialImage:
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"{name} must be a NIfTI image")
    result = normalize_image_to_ras(image, source_name=name)
    if len(result.shape) != 3:
        raise ValueError(f"{name} must be a 3D scalar image")
    return result


def _make_output(
    data: np.ndarray,
    reference: nib.spatialimages.SpatialImage,
) -> nib.Nifti1Image:
    header = reference.header.copy()
    header.set_data_dtype(data.dtype)
    header.set_slope_inter(None, None)
    return nib.Nifti1Image(data, reference.affine, header=header)
