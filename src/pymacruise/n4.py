"""ANTsPy N4 correction with the selected V3 command-line defaults."""

from __future__ import annotations

from dataclasses import dataclass
import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras


class N4DependencyError(RuntimeError):
    """Raised when the optional ANTsPy dependency is not usable."""


@dataclass(frozen=True)
class N4Config:
    """N4 settings matching the selected wrapper's unadorned CLI invocation.

    The MATLAB wrapper calls ``N4BiasFieldCorrection -i input -o output``.
    The ANTs command-line defaults are 50 iterations at four levels, zero
    convergence threshold, shrink factor 4, and one B-spline mesh element per
    axis. ANTsPy's tolerance default differs, so it is set explicitly here.
    """

    shrink_factor: int = 4
    iterations: tuple[int, ...] = (50, 50, 50, 50)
    convergence_tolerance: float = 0.0
    spline_mesh: tuple[int, int, int] = (1, 1, 1)
    rescale_intensities: bool = False

    def __post_init__(self) -> None:
        if self.shrink_factor < 1:
            raise ValueError("N4 shrink factor must be a positive integer")
        if not self.iterations or any(count < 1 for count in self.iterations):
            raise ValueError("N4 requires one or more positive iteration counts")
        if not np.isfinite(self.convergence_tolerance) or self.convergence_tolerance < 0:
            raise ValueError("N4 convergence tolerance must be finite and non-negative")
        if len(self.spline_mesh) != 3 or any(size < 1 for size in self.spline_mesh):
            raise ValueError("N4 spline mesh must contain three positive integer sizes")


def _nib_to_ants(image: nib.spatialimages.SpatialImage):
    """Create an ANTsPy image while explicitly converting RAS geometry to LPS."""
    try:
        import ants
    except (ImportError, OSError) as exc:
        raise N4DependencyError(
            "ANTsPy is required for N4 correction; install PyMaCRUISE with its 'ants' extra"
        ) from exc

    affine = np.asarray(image.affine, dtype=np.float64)
    spacing = nib.affines.voxel_sizes(affine).astype(np.float64)
    ras_to_lps = np.diag([-1.0, -1.0, 1.0])
    direction_lps = ras_to_lps @ (affine[:3, :3] / spacing[np.newaxis, :])
    origin_lps = ras_to_lps @ affine[:3, 3]
    values = np.asarray(image.dataobj, dtype=np.float32)
    ants_image = ants.from_numpy(
        values,
        origin=tuple(float(v) for v in origin_lps),
        spacing=tuple(float(v) for v in spacing),
        direction=direction_lps,
    )
    return ants, ants_image


def n4_bias_correct(
    image: nib.spatialimages.SpatialImage,
    *,
    config: N4Config = N4Config(),
) -> nib.spatialimages.SpatialImage:
    """Correct one 3-D intensity image and return it on its unchanged RAS grid.

    The wrapper mirrors the source command's omission of an explicit mask and
    uses the ANTs CLI-equivalent N4 defaults. Callers needing a mask or a finer
    spline mesh must supply a separately validated source-compatible profile.
    """
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("N4 input must be a NIfTI image")
    if len(image.shape) != 3:
        raise ValueError(f"N4 input must be a 3D scalar volume; got {image.shape}")
    normalized = normalize_image_to_ras(image, source_name="N4 input")
    values = np.asarray(normalized.dataobj, dtype=np.float32)
    if not np.isfinite(values).all():
        raise ValueError("N4 input must contain finite voxel values")
    if np.any(values < 0):
        raise ValueError("N4 input contains negative values; the ANTs log-domain method requires non-negative intensities")
    if not np.any(values > 0):
        raise ValueError("N4 input must contain at least one positive intensity")

    ants, ants_image = _nib_to_ants(normalized)
    try:
        corrected = ants.n4_bias_field_correction(
            ants_image,
            mask=None,
            rescale_intensities=config.rescale_intensities,
            shrink_factor=config.shrink_factor,
            convergence={
                "iters": list(config.iterations),
                "tol": float(config.convergence_tolerance),
            },
            spline_param=list(config.spline_mesh),
            return_bias_field=False,
            verbose=False,
        )
    except Exception as exc:
        raise RuntimeError(f"ANTsPy N4 bias correction failed: {exc}") from exc

    corrected_data = np.asarray(corrected.numpy(), dtype=np.float32)
    if corrected_data.shape != normalized.shape or not np.isfinite(corrected_data).all():
        raise RuntimeError("ANTsPy N4 returned invalid data or changed the input grid")
    header = normalized.header.copy()
    header.set_data_dtype(np.float32)
    header.set_slope_inter(None, None)
    result = nib.Nifti1Image(corrected_data, normalized.affine, header=header)
    qform, qcode = normalized.get_qform(coded=True)
    _, scode = normalized.get_sform(coded=True)
    result.set_qform(qform if qcode else None, code=int(qcode), strip_shears=True)
    result.set_sform(normalized.affine, code=int(scode) if scode else 1)
    return result
