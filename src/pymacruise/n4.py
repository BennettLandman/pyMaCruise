"""ANTsPy N4 correction with explicit modern and historical mask semantics.

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
    """Shared N4 schedule with an explicit ANTs-generation mask policy.

    The MATLAB wrapper calls ``N4BiasFieldCorrection -i input -o output``.
    The ANTs command-line defaults are 50 iterations at four levels, zero
    convergence threshold, shrink factor 4, one B-spline mesh element per
    axis, and intensity rescaling disabled. ANTs CLI only enables its
    rescaling post-step when an explicit mask is supplied; the wrapper supplies
    none. Modern N4 uses the entire image in that case. Historical ANTs 2.1.0
    instead estimates an automatic 200-bin Otsu mask. The default follows the
    current unmasked backend; BASE can select ``legacy-otsu`` explicitly.
    ANTsPy's tolerance default differs, so it is set explicitly.
    """

    shrink_factor: int = 4
    iterations: tuple[int, ...] = (50, 50, 50, 50)
    convergence_tolerance: float = 0.0
    spline_mesh: tuple[int, int, int] = (1, 1, 1)
    rescale_intensities: bool = False
    estimation_mask: str = "all"

    def __post_init__(self) -> None:
        if self.estimation_mask not in ("all", "legacy-otsu"):
            raise ValueError("N4 estimation_mask must be all or legacy-otsu")
        if self.estimation_mask == "legacy-otsu" and self.rescale_intensities:
            raise ValueError("Historical automatic Otsu N4 does not rescale intensities")
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
        ants, converted = nifti_to_ants(image)
        # ITK's NIfTI reader takes spacing from pixdim, while the shared
        # geometry adapter derives it from affine column lengths. Their tiny
        # floating-point differences can change N4's iterative spline fit.
        # Preserve the exact stored spacing used by a bare ANTs CLI call.
        converted.set_spacing(tuple(float(value) for value in image.header.get_zooms()[:3]))
        return ants, converted
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

    The ``all`` policy follows the current bare ANTs backend. ``legacy-otsu``
    reproduces historical automatic mask estimation while correcting all
    voxels, as the older unmasked command did. Both policies retain the input
    acquisition lattice; working-grid conversion belongs after estimation.
    """
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("N4 input must be a NIfTI image")
    if len(image.shape) != 3:
        raise ValueError(f"N4 input must be a 3D scalar volume; got {image.shape}")
    normalized = normalize_image_to_ras(image, source_name="N4 input")
    values = np.asarray(normalized.dataobj, dtype=np.float32)
    if not np.isfinite(values).all():
        raise ValueError("N4 input must contain finite voxel values")
    if not np.any(values > 0):
        raise ValueError("N4 input must contain at least one positive intensity")

    # The source wrapper invokes the ANTs CLI without specifying a mask or
    # algorithm parameters. These values reproduce that wrapper's CLI defaults
    # while overriding ANTsPy's different tolerance default; they are workflow
    # settings, not numerical values prescribed by the N4 paper. The N4
    # histogram/spline estimation itself is delegated to ANTsPy.
    ants, ants_image = _nib_to_ants(normalized)
    estimation_mask = None
    legacy_otsu = config.estimation_mask == "legacy-otsu"
    if legacy_otsu:
        from pymacruise.n4_legacy_otsu import legacy_n4_otsu_mask

        mask_values, otsu_threshold = legacy_n4_otsu_mask(values)
        estimation_mask = ants_image.new_image_like(mask_values.astype(np.float32))
    try:
        # Invoke the published N4 bias-estimation implementation through the
        # ANTsPy API: the paper's iterative intensity correction and smooth
        # B-spline bias-field model are in ANTs, not duplicated here.
        corrected = ants.n4_bias_field_correction(
            ants_image,
            mask=estimation_mask,
            rescale_intensities=config.rescale_intensities,
            shrink_factor=config.shrink_factor,
            convergence={
                "iters": list(config.iterations),
                "tol": float(config.convergence_tolerance),
            },
            spline_param=list(config.spline_mesh),
            return_bias_field=legacy_otsu,
            verbose=False,
        )
    except Exception as exc:
        raise RuntimeError(f"ANTsPy N4 bias correction failed: {exc}") from exc

    corrected_data = np.asarray(corrected.numpy(), dtype=np.float32)
    if legacy_otsu:
        # Historical N4's automatically generated mask controls estimation;
        # its final correction still covers the whole image. Modern ANTs
        # restores input values outside an explicitly supplied mask. Apply
        # the estimated bias field ourselves to preserve the historical rule.
        if corrected_data.shape != values.shape or not np.isfinite(corrected_data).all() or np.any(corrected_data <= 0):
            raise RuntimeError("ANTsPy N4 returned an invalid bias field")
        corrected_data = np.divide(values, corrected_data, dtype=np.float32)
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
    if legacy_otsu:
        result.extra["n4_estimation_mask"] = {
            "strategy": "legacy-otsu", "histogram_bins": 200,
            "threshold": otsu_threshold, "estimation_voxels": int(mask_values.sum()),
            "output_correction": "entire-image",
        }
    return result
