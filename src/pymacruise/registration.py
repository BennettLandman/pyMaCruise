"""ANTsPy registration and physical-grid transfer for MaCRUISE.

The selected MATLAB adapters call NiftyReg ``reg_aladin`` with the template
as reference and the subject as floating image. T1/FLAIR and MNI alignment
use rigid-only mode in the selected runners. The NiftyReg method is based on
symmetric block matching (Modat et al., 2014,
https://doi.org/10.1117/1.JMI.1.2.024003). ``reg_resample`` maps floating
images into the reference grid; labels explicitly request nearest-neighbor
interpolation, while unspecified continuous images use NiftyReg's cubic
spline default. This module instead calls ANTsPy/ITK rigid or affine
registration and transform application (ANTsPy API:
https://antspy.readthedocs.io/en/latest/api/ants.registration.html). Its
metric, optimizer, transform implementation, and spline resampling are a
software substitution; they are not claimed to reproduce NiftyReg.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, replace
import os
from pathlib import Path
from typing import Literal, Sequence

import nibabel as nib
import numpy as np

from pymacruise.ants_utils import ANTsDependencyError, nifti_to_ants
from pymacruise.image_io import normalize_image_to_ras
from pymacruise.surfaces import SurfaceMesh, _signed_mesh_volume


class ANTsRegistrationError(RuntimeError):
    """Raised when ANTsPy registration or resampling fails."""


@dataclass(frozen=True)
class RegistrationConfig:
    """Explicit ANTsPy settings for global rigid/affine alignment.

    The ANTsPy API exposes the metric, sampling, iterations, and multiresolution
    schedule used below; the selected values are this adapter's configuration,
    not settings prescribed by the ANTs metric-evaluation paper or a claim of
    matching NiftyReg Aladin (Modat et al., 2014,
    https://doi.org/10.1117/1.JMI.1.2.024003).
    """

    transform_type: Literal["Rigid", "Affine"] = "Rigid"
    metric: Literal["mattes", "meansquares", "GC"] = "mattes"
    metric_bins: int = 32
    sampling_rate: float = 0.2
    iterations: tuple[int, ...] = (1000, 500, 250, 0)
    shrink_factors: tuple[int, ...] = (6, 4, 2, 1)
    smoothing_sigmas: tuple[float, ...] = (3.0, 2.0, 1.0, 0.0)
    deterministic: bool = True
    random_seed: int = 123
    single_precision: bool = True

    def __post_init__(self) -> None:
        if self.transform_type not in {"Rigid", "Affine"}:
            raise ValueError("Registration transform_type must be Rigid or Affine")
        if self.metric_bins < 2:
            raise ValueError("Registration metric_bins must be at least 2")
        if not 0 < self.sampling_rate <= 1:
            raise ValueError("Registration sampling_rate must be in (0, 1]")
        if not self.iterations or len(self.iterations) != len(self.shrink_factors):
            raise ValueError("Registration iteration and shrink schedules must have equal length")
        if len(self.iterations) != len(self.smoothing_sigmas):
            raise ValueError("Registration iteration and smoothing schedules must have equal length")
        if any(value < 0 for value in self.iterations):
            raise ValueError("Registration iterations must be nonnegative")
        if any(value < 1 for value in self.shrink_factors):
            raise ValueError("Registration shrink factors must be positive")
        if any(not np.isfinite(value) or value < 0 for value in self.smoothing_sigmas):
            raise ValueError("Registration smoothing sigmas must be finite and nonnegative")


@dataclass(frozen=True)
class RegistrationResult:
    """ANTs image-resampling transforms, plus a moving image on the fixed grid."""

    forward_transforms: tuple[str, ...]
    inverse_transforms: tuple[str, ...]
    warped_moving: nib.Nifti1Image
    transform_type: str


def _validated_ras_image(
    image: nib.spatialimages.SpatialImage,
    name: str,
) -> nib.spatialimages.SpatialImage:
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"{name} must be a NIfTI image")
    if len(image.shape) != 3:
        raise ValueError(f"{name} must be a 3D scalar image")
    result = normalize_image_to_ras(image, source_name=name)
    values = np.asarray(result.dataobj)
    if not np.isfinite(values).all():
        raise ValueError(f"{name} contains NaN or infinite values")
    return result


def _nifti_on_fixed_grid(
    ants_image,
    fixed: nib.spatialimages.SpatialImage,
    *,
    dtype: np.dtype | str,
) -> nib.Nifti1Image:
    values = np.asarray(ants_image.numpy(), dtype=np.dtype(dtype))
    if values.shape != fixed.shape or not np.isfinite(values).all():
        raise ANTsRegistrationError("ANTsPy output is non-finite or changed the fixed grid")
    header = fixed.header.copy()
    header.set_data_dtype(values.dtype)
    header.set_slope_inter(None, None)
    result = nib.Nifti1Image(values, fixed.affine, header=header)
    qform, qcode = fixed.get_qform(coded=True)
    _, scode = fixed.get_sform(coded=True)
    result.set_qform(qform if qcode else None, code=int(qcode), strip_shears=True)
    result.set_sform(fixed.affine, code=int(scode) if scode else 1)
    return result


def _validate_transform_list(transform_list: Sequence[str | Path]) -> tuple[str, ...]:
    if not transform_list:
        raise ValueError("At least one ANTs transform is required")
    transforms = tuple(str(Path(path)) for path in transform_list)
    missing = [path for path in transforms if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError(f"ANTs transform file does not exist: {missing[0]}")
    return transforms


@contextmanager
def _deterministic_ants(ants, enabled: bool, seed: int):
    if not enabled:
        yield
        return
    previous_deterministic = ants.config._deterministic
    previous_seed = ants.config._random_seed
    previous_threads = os.environ.get("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS")
    previous_ants_seed = os.environ.get("ANTS_RANDOM_SEED")
    ants.config.set_ants_deterministic(True, seed_value=seed)
    if seed is not None:
        # ANTs registration reads its C++ random seed from this environment
        # variable. ANTsPy's NumPy/Python seeds alone do not control the
        # random metric sampling used by Rigid/Affine registration.
        os.environ["ANTS_RANDOM_SEED"] = str(seed)
    try:
        yield
    finally:
        ants.config.set_ants_deterministic(previous_deterministic, seed_value=previous_seed)
        if previous_ants_seed is None:
            os.environ.pop("ANTS_RANDOM_SEED", None)
        else:
            os.environ["ANTS_RANDOM_SEED"] = previous_ants_seed
        if previous_threads is None:
            os.environ.pop("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS", None)
        else:
            os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = previous_threads


@contextmanager
def _ants_runtime_environment(enabled: bool, seed: int):
    """Set seed and thread count before ANTsPy/ITK initializes its runtime."""
    if not enabled or seed is None:
        yield
        return
    previous_seed = os.environ.get("ANTS_RANDOM_SEED")
    previous_threads = os.environ.get("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS")
    os.environ["ANTS_RANDOM_SEED"] = str(seed)
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
    try:
        yield
    finally:
        if previous_seed is None:
            os.environ.pop("ANTS_RANDOM_SEED", None)
        else:
            os.environ["ANTS_RANDOM_SEED"] = previous_seed
        if previous_threads is None:
            os.environ.pop("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS", None)
        else:
            os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = previous_threads


def register_to_reference(
    fixed: nib.spatialimages.SpatialImage,
    moving: nib.spatialimages.SpatialImage,
    *,
    output_prefix: str | Path,
    fixed_mask: nib.spatialimages.SpatialImage | None = None,
    config: RegistrationConfig = RegistrationConfig(),
) -> RegistrationResult:
    # ANTs/ITK can read ANTS_RANDOM_SEED while the native runtime is loading,
    # before ANTsPy's config seed is applied. Set it before converting images
    # to ANTs objects so repeated process-level registrations use the same seed.
    with _ants_runtime_environment(config.deterministic, config.random_seed):
        return _register_to_reference(
            fixed,
            moving,
            output_prefix=output_prefix,
            fixed_mask=fixed_mask,
            config=config,
        )


def _register_to_reference(
    fixed: nib.spatialimages.SpatialImage,
    moving: nib.spatialimages.SpatialImage,
    *,
    output_prefix: str | Path,
    fixed_mask: nib.spatialimages.SpatialImage | None = None,
    config: RegistrationConfig = RegistrationConfig(),
) -> RegistrationResult:
    """Estimate a moving-to-fixed rigid/affine transform with ANTsPy.

    ``fixed`` is the legacy ``-ref`` image; ``moving`` is ``-flo``. The
    returned forward transform list follows ANTs image-resampling convention;
    point application uses that list in the opposite direction.
    Output prefix is required so transform files have a caller-controlled,
    persistent location. The optional reference-space mask applies to all
    registration stages, using the same fixed/reference-space mask role as
    NiftyReg ``-rmask``. The backend behavior remains ANTsPy, not Aladin's
    symmetric block-matching method (Modat et al., 2014,
    https://doi.org/10.1117/1.JMI.1.2.024003); the ANTs API documents fixed,
    moving, mask, metric, and schedule arguments at
    https://antspy.readthedocs.io/en/latest/api/ants.registration.html.
    """
    fixed_ras = _validated_ras_image(fixed, "fixed registration image")
    moving_ras = _validated_ras_image(moving, "moving registration image")
    fixed_ants_pair = nifti_to_ants(fixed_ras)
    moving_ants_pair = nifti_to_ants(moving_ras)
    ants, fixed_ants = fixed_ants_pair
    if ants is not moving_ants_pair[0]:
        raise RuntimeError("ANTsPy module changed during registration setup")
    moving_ants = moving_ants_pair[1]

    mask_ants = None
    if fixed_mask is not None:
        mask_ras = _validated_ras_image(fixed_mask, "fixed registration mask")
        if mask_ras.shape != fixed_ras.shape or not np.allclose(
            mask_ras.affine, fixed_ras.affine, rtol=0, atol=1e-4
        ):
            raise ValueError("Fixed registration mask must share the fixed RAS grid")
        mask_data = (np.asarray(mask_ras.dataobj) != 0).astype(np.uint8)
        mask_ants = ants.from_numpy(
            mask_data,
            origin=fixed_ants.origin,
            spacing=fixed_ants.spacing,
            direction=fixed_ants.direction,
        )
        if not np.any(mask_data):
            raise ValueError("Fixed registration mask is empty")

    prefix = Path(output_prefix)
    if not prefix.name:
        raise ValueError("Registration output_prefix must have a filename prefix")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    try:
        with _deterministic_ants(ants, config.deterministic, config.random_seed):
            outputs = ants.registration(
                fixed=fixed_ants,
                moving=moving_ants,
                type_of_transform=config.transform_type,
                outprefix=str(prefix),
                mask=mask_ants,
                mask_all_stages=mask_ants is not None,
                aff_metric=config.metric,
                aff_sampling=config.metric_bins,
                aff_random_sampling_rate=config.sampling_rate,
                aff_iterations=config.iterations,
                aff_shrink_factors=config.shrink_factors,
                aff_smoothing_sigmas=config.smoothing_sigmas,
                singleprecision=config.single_precision,
                verbose=False,
            )
    except Exception as exc:
        raise ANTsRegistrationError(f"ANTsPy {config.transform_type} registration failed: {exc}") from exc

    forward = _validate_transform_list(outputs.get("fwdtransforms", ()))
    inverse = _validate_transform_list(outputs.get("invtransforms", ()))
    warped = apply_image_transform(
        fixed_ras,
        moving_ras,
        forward,
        interpolation="intensity",
    )
    return RegistrationResult(forward, inverse, warped, config.transform_type)


def register_flair_to_t1(
    t1: nib.spatialimages.SpatialImage,
    flair: nib.spatialimages.SpatialImage,
    *,
    output_prefix: str | Path,
    fixed_mask: nib.spatialimages.SpatialImage | None = None,
    config: RegistrationConfig = RegistrationConfig(),
) -> RegistrationResult:
    """Register FLAIR (moving) to the T1 reference with a rigid transform.

    This is the selected MaCRUISEplus `FlairRegistration.m` contract. ANTsPy
    settings remain caller-configurable, while transform type is pinned to
    rigid-only as in the source wrapper.
    """
    if config.transform_type != "Rigid":
        raise ValueError("MaCRUISEplus FLAIR-to-T1 registration must be rigid-only")
    return register_to_reference(
        t1,
        flair,
        output_prefix=output_prefix,
        fixed_mask=fixed_mask,
        config=config,
    )


def apply_image_transform(
    fixed: nib.spatialimages.SpatialImage,
    moving: nib.spatialimages.SpatialImage,
    transform_list: Sequence[str | Path],
    *,
    interpolation: Literal["intensity", "levelset", "label"] = "intensity",
    default_value: float = 0.0,
) -> nib.Nifti1Image:
    """Resample a moving NIfTI into the fixed RAS grid using ANTs transforms.

    Intensity and level-set values use ANTsPy's ``bSpline`` interpolator;
    labels use ``nearestNeighbor``. The selected NiftyReg wrappers use
    cubic-spline default interpolation for continuous images and ``-inter 0``
    for labels. The choices follow the same interpolation categories, but the
    backend kernels, boundary semantics, precision, and outputs are not
    presumed identical. ANTsPy API: https://antspyx.readthedocs.io/en/stable/registration.html;
    NiftyReg default: https://github.com/SuperElastix/niftyreg/blob/master/reg-apps/reg_resample.cpp.
    """
    if interpolation not in {"intensity", "levelset", "label"}:
        raise ValueError("Interpolation must be intensity, levelset, or label")
    if not np.isfinite(default_value):
        raise ValueError("Outside-grid default value must be finite")
    fixed_ras = _validated_ras_image(fixed, "fixed image")
    moving_ras = _validated_ras_image(moving, "moving image")
    transforms = _validate_transform_list(transform_list)
    ants, fixed_ants = nifti_to_ants(fixed_ras)
    _, moving_ants = nifti_to_ants(moving_ras)
    if interpolation == "label":
        values = np.asarray(moving_ras.dataobj)
        if not np.equal(values, np.rint(values)).all():
            raise ValueError("Label resampling requires integer-valued moving data")
        ants_moving = ants.from_numpy(
            values.astype(np.float32),
            origin=moving_ants.origin,
            spacing=moving_ants.spacing,
            direction=moving_ants.direction,
        )
        # Preserve the wrappers' nearest-neighbor label-resampling intent.
        # This is a corresponding ANTsPy interpolator choice, not a claim that
        # every NiftyReg/ITK boundary or tie case has identical output.
        interpolator = "nearestNeighbor"
        output_dtype = moving_ras.get_data_dtype()
    else:
        ants_moving = moving_ants
        # reg_resample defaults to cubic spline, and ANTsPy offers bSpline.
        # Keep this as an explicit backend substitution: names do not establish
        # identical spline kernels or boundary behavior.
        interpolator = "bSpline"
        output_dtype = np.float32
    try:
        warped = ants.apply_transforms(
            fixed=fixed_ants,
            moving=ants_moving,
            transformlist=list(transforms),
            interpolator=interpolator,
            defaultvalue=float(default_value),
            singleprecision=True,
            verbose=False,
        )
    except Exception as exc:
        raise ANTsRegistrationError(f"ANTsPy transform application failed: {exc}") from exc
    return _nifti_on_fixed_grid(warped, fixed_ras, dtype=output_dtype)


def apply_transforms_to_points_ras(
    points_ras: np.ndarray,
    transform_list: Sequence[str | Path],
) -> np.ndarray:
    """Transform moving-space RAS+ millimetre points into the target RAS frame.

    ANTs transform lists used for image resampling map fixed-grid sample points
    back into the moving image. Point mapping runs in the opposite direction,
    so each transform is inverted while evaluating the composite. This follows
    ANTsPy/ITK transform-list conventions (see
    https://antspyx.readthedocs.io/en/stable/registration.html); RAS/LPS
    conversion and affine-only point handling are this project's adapter logic.
    """
    points = np.asarray(points_ras, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("Points must be a finite (N, 3) RAS coordinate array")
    if len(points) == 0:
        return points.copy()
    transforms = _validate_transform_list(transform_list)
    try:
        import ants
    except (ImportError, OSError) as exc:
        raise ANTsDependencyError(
            "ANTsPy is required for point transforms; install PyMaCRUISE with its 'ants' extra"
        ) from exc

    if any(Path(path).suffix != ".mat" for path in transforms):
        raise ValueError("RAS point transforms currently support affine ANTs .mat files only")
    # The C++ antsApplyTransformsToPoints entry point has crashed with otherwise
    # valid point tables in some ANTsPy builds. Read affine transforms and use
    # ANTsPy's per-point API. Invert the image-resampling maps and reverse list
    # order to match ITK composite transform stack evaluation.
    # ANTs/ITK transforms use LPS physical coordinates; PyMaCRUISE exposes RAS.
    points_lps = points.copy()
    points_lps[:, :2] *= -1
    try:
        txs = [
            ants.invert_ants_transform(ants.read_transform(path, precision="double"))
            for path in reversed(transforms)
        ]
        output_lps = np.asarray(
            [
                _apply_point_transform_chain(ants, txs, point)
                for point in points_lps
            ],
            dtype=np.float64,
        )
    except Exception as exc:
        raise ANTsRegistrationError(f"ANTsPy point transform failed: {exc}") from exc
    output_lps[:, :2] *= -1
    if not np.isfinite(output_lps).all():
        raise ANTsRegistrationError("ANTsPy point transform returned non-finite coordinates")
    return output_lps


def _apply_point_transform_chain(ants, transforms, point: np.ndarray) -> tuple[float, ...]:
    result = tuple(float(value) for value in point)
    for transform in transforms:
        result = ants.apply_ants_transform_to_point(transform, result)
    return result


def apply_transforms_to_surface_ras(
    mesh: SurfaceMesh,
    transform_list: Sequence[str | Path],
) -> SurfaceMesh:
    """Transform a scientific mesh's vertices while retaining faces/data."""
    vertices = apply_transforms_to_points_ras(mesh.vertices, transform_list)
    faces = mesh.faces.copy()
    before_volume = _signed_mesh_volume(mesh.vertices, faces)
    after_volume = _signed_mesh_volume(vertices, faces)
    if before_volume != 0 and after_volume != 0 and np.sign(before_volume) != np.sign(after_volume):
        faces = faces[:, ::-1]
    return replace(mesh, vertices=vertices, faces=faces, gifti_dataspace=0)
