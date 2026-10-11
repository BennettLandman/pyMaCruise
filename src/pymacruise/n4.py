"""ANTsPy N4 correction with the selected V3 command-line defaults.

The N4 estimation algorithm is provided by ANTsPy/ANTs; this module adapts
image geometry, wrapper-selected parameters, and output metadata rather than
reimplementing the method (Tustison et al., 2010; see
``docs/papers/n4.md``).
"""

from __future__ import annotations

from dataclasses import dataclass
import nibabel as nib
import numpy as np

from pymacruise.ants_utils import ANTsDependencyError, nifti_to_ants
from pymacruise.image_io import normalize_image_to_ras


class N4DependencyError(ANTsDependencyError):
    """Raised when the optional ANTsPy dependency is not usable."""


@dataclass(frozen=True)
class N4Config:
    """N4 settings matching the selected wrapper's unadorned CLI invocation.

    The MATLAB wrapper calls ``N4BiasFieldCorrection -i input -o output``.
    The ANTs command-line defaults are 50 iterations at four levels, zero
    convergence threshold, shrink factor 4, one B-spline mesh element per
    axis, and intensity rescaling enabled. ANTsPy's tolerance and rescaling
    defaults differ, so both are set explicitly here.
    """

    shrink_factor: int = 4
    iterations: tuple[int, ...] = (50, 50, 50, 50)
    convergence_tolerance: float = 0.0
    spline_mesh: tuple[int, int, int] = (1, 1, 1)
    rescale_intensities: bool = True

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
        return nifti_to_ants(image)
    except ANTsDependencyError as exc:
        raise N4DependencyError(
            "ANTsPy is required for N4 correction; install PyMaCRUISE with its 'ants' extra"
        ) from exc


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

    # The source wrapper invokes the ANTs CLI without specifying a mask or
    # algorithm parameters. These values reproduce that wrapper's CLI defaults
    # while overriding ANTsPy's different tolerance default; they are workflow
    # settings, not numerical values prescribed by the N4 paper. The N4
    # histogram/spline estimation itself is delegated to ANTsPy.
    ants, ants_image = _nib_to_ants(normalized)
    try:
        # Invoke the published N4 bias-estimation implementation through the
        # ANTsPy API: the paper's iterative intensity correction and smooth
        # B-spline bias-field model are in ANTs, not duplicated here.
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
