from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_approximate_centroids import update_default_toads_approximate_centroids
from pymacruise.toads_approximate_variances import compute_default_toads_approximate_variances
from pymacruise.toads_default_approximate_iteration import run_default_toads_approximate_iteration
from pymacruise.toads_default_approximate_pass import toads_default_approximate_membership_pass


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_approximate_pass"
VARIANCE_ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_approximate_variances"
CENTROID_ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_approximate_centroids"


def _nifti(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj).astype(np.float32, copy=True)


def _txt(name: str) -> np.ndarray:
    return np.loadtxt(ROOT / f"{name}.txt", dtype=np.float32).reshape(-1)


def _inputs() -> dict[str, np.ndarray | float | None]:
    params = _txt("parameters")
    return {
        "images": _nifti("images"),
        "memberships": _nifti("memberships"),
        "best_class_indices": _nifti("best").astype(np.uint8),
        "segmentation": np.asarray(
            nib.load(str(VARIANCE_ROOT / "segmentation.nii.gz")).dataobj
        ).astype(np.uint8),
        "centroids": _txt("centroids").reshape(1, 11),
        "intensity_field": _nifti("field"),
        "optimized_exponents": np.loadtxt(
            CENTROID_ROOT / "optimized_exponents.txt", dtype=np.float32
        ).reshape(-1, 3)[:1],
        "modality_weights": _txt("modality"),
        "exponent_by_class": _txt("exponents"),
        "class_variance_weights": _txt("variance_weights").reshape(1, 11),
        "transformed_shapes": _nifti("transformed_shapes"),
        "direct_relation": _nifti("direct_relation") != 0,
        "smoothing_priors": _txt("smoothing_priors"),
        "atlas_shape_priors": _txt("shape_priors"),
        "prior_scale": float(params[0]),
        "prior_factor": float(params[1]),
        "prior_coefficient": float(params[2]),
        "smoothing": float(params[3]),
    }


def test_approximate_iteration_composes_centroids_variances_and_memberships() -> None:
    values = _inputs()
    image = values["images"]
    memberships = values["memberships"]
    best = values["best_class_indices"]
    segmentation = values["segmentation"]
    centers = values["centroids"]
    field = values["intensity_field"]
    exponents = values["optimized_exponents"]
    snapshots = {
        key: value.copy()
        for key, value in values.items()
        if isinstance(value, np.ndarray)
    }

    result = run_default_toads_approximate_iteration(**values)

    centroid = update_default_toads_approximate_centroids(
        image, memberships, best, centers, field, exponents,
    )
    variance = compute_default_toads_approximate_variances(
        image, memberships, best, segmentation, centroid.centroids,
        intensity_field=field,
    )
    expected_membership = toads_default_approximate_membership_pass(
        image,
        memberships,
        best,
        centroid.centroids,
        variance.standard_deviations,
        values["class_variance_weights"],
        values["modality_weights"],
        values["exponent_by_class"],
        values["transformed_shapes"],
        values["direct_relation"],
        values["smoothing_priors"],
        values["atlas_shape_priors"],
        prior_scale=values["prior_scale"],
        prior_factor=values["prior_factor"],
        prior_coefficient=values["prior_coefficient"],
        smoothing=values["smoothing"],
        intensity_field=field,
    )
    np.testing.assert_array_equal(result.centroids.view(np.uint32), centroid.centroids.view(np.uint32))
    np.testing.assert_array_equal(result.centroid_numerators.view(np.uint32), centroid.numerators.view(np.uint32))
    np.testing.assert_array_equal(result.centroid_denominators.view(np.uint32), centroid.denominators.view(np.uint32))
    np.testing.assert_array_equal(
        result.standard_deviations.view(np.uint32), variance.standard_deviations.view(np.uint32)
    )
    np.testing.assert_array_equal(
        result.variance_numerators.view(np.uint32), variance.numerators.view(np.uint32)
    )
    np.testing.assert_array_equal(
        result.variance_denominators.view(np.uint32), variance.denominators.view(np.uint32)
    )
    np.testing.assert_array_equal(
        result.membership.memberships.view(np.uint32), expected_membership.memberships.view(np.uint32)
    )
    np.testing.assert_array_equal(
        result.membership.best_class_indices, expected_membership.best_class_indices
    )
    np.testing.assert_array_equal(
        result.membership.denominator.view(np.uint32), expected_membership.denominator.view(np.uint32)
    )
    assert result.membership.objective.view(np.uint32) == expected_membership.objective.view(np.uint32)
    for key, before in snapshots.items():
        np.testing.assert_array_equal(values[key], before)
