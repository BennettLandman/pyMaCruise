from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_approximate_pass import (
    toads_default_approximate_membership_pass,
)

ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_approximate_pass"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsDefaultApproximateMembershipPassHarness.java"
HARNESS_SHA256 = "df02e10efc85007255cb00cdbb692b44197a562e1f2bf3c78c8db65302086007"
SOURCE_HASHES = {
    "segmentation_sha256": "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884",
    "algorithm_sha256": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
    "wrapper_sha256": "d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115",
    "descriptor_sha256": "d43efa5a002b0c1a838449ba7ca1d3bf35fab5e3b575474a71aadfe475dcc257",
    "numerics_sha256": "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb",
    "powertable_sha256": "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678",
}


def _nifti(name: str) -> np.ndarray:
    volume = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
    np.testing.assert_allclose(
        nib.affines.voxel_sizes(volume.affine), (.85, 1.15, 1.55), rtol=0, atol=1e-6
    )
    return np.asarray(volume.dataobj).astype(np.float32, copy=True)


def _txt(name: str) -> np.ndarray:
    return np.loadtxt(ROOT / f"{name}.txt", dtype=np.float32).reshape(-1)


def _bits(values: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32).reshape(shape)


def _run(*, use_field: bool = True, overrides: dict | None = None):
    params = np.loadtxt(ROOT / "parameters.txt", dtype=np.float32)
    values = dict(
        image_values=_nifti("images"),
        initial_memberships=_nifti("memberships"),
        best_class_indices=_nifti("best").astype(np.uint8),
        centroids=_txt("centroids").reshape(1, 11),
        standard_deviations=_txt("standard_deviations").reshape(1, 11),
        class_variance_weights=_txt("variance_weights").reshape(1, 11),
        modality_weights=_txt("modality"),
        exponent_by_class=_txt("exponents"),
        transformed_shapes=_nifti("transformed_shapes"),
        direct_relation=_nifti("direct_relation") != 0,
        smoothing_priors=_txt("smoothing_priors"),
        atlas_shape_priors=_txt("shape_priors"),
        prior_scale=float(params[0]),
        prior_factor=float(params[1]),
        prior_coefficient=float(params[2]),
        smoothing=float(params[3]),
        intensity_field=_nifti("field") if use_field else None,
    )
    if overrides:
        values.update(overrides)
    return toads_default_approximate_membership_pass(**values)


def test_alternate_approximate_membership_pass_matches_java_bits() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert all(JAVA[key] == value for key, value in SOURCE_HASHES.items())
    assert JAVA["validation_profile"] == (
        "alternate approx=4 constructor path; MedicAlgorithmToads08 shorter constructor forces useApprox=false"
    )
    shape = tuple(JAVA["shape"][:3])
    assert JAVA["shape"][3] == 1 and JAVA["classes"] == 11 and JAVA["approx"] == 4
    result = _run(use_field=True)
    np.testing.assert_array_equal(result.image_variances.view(np.int32), _bits(JAVA["field_image_variance_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(result.image_variance_numerator.view(np.int32), _bits(JAVA["field_image_variance_num_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(result.image_variance_denominator.view(np.int32), _bits(JAVA["field_image_variance_den_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(result.unnormalized_image_variance.view(np.int32), _bits(JAVA["field_image_variance_raw_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(result.image_standard_deviation.view(np.int32), _bits(JAVA["field_image_std_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(result.centroid_similarity.view(np.int32), _bits(JAVA["field_similarity_bits"], (11, 11)).view(np.int32))
    np.testing.assert_array_equal(
        result.memberships.transpose(2, 1, 0, 3).reshape(-1).view(np.int32),
        _bits(JAVA["field_membership_bits"], (np.prod(shape) * 4,)).view(np.int32),
    )
    np.testing.assert_array_equal(
        result.best_class_indices.transpose(2, 1, 0, 3).reshape(-1),
        np.asarray(JAVA["field_best"], dtype=np.uint8),
    )
    np.testing.assert_array_equal(result.denominator.view(np.int32), _bits(JAVA["field_denominator_bits"], shape).view(np.int32))
    assert result.objective.view(np.uint32) == np.asarray(JAVA["field_objective_bits"], dtype=np.uint32)

    raw = _run(use_field=False)
    np.testing.assert_array_equal(raw.image_variances.view(np.int32), _bits(JAVA["raw_image_variance_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(raw.image_variance_numerator.view(np.int32), _bits(JAVA["raw_image_variance_num_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(raw.image_variance_denominator.view(np.int32), _bits(JAVA["raw_image_variance_den_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(raw.unnormalized_image_variance.view(np.int32), _bits(JAVA["raw_image_variance_raw_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(raw.image_standard_deviation.view(np.int32), _bits(JAVA["raw_image_std_bits"], (1,)).view(np.int32))
    np.testing.assert_array_equal(
        raw.memberships.transpose(2, 1, 0, 3).reshape(-1).view(np.int32),
        _bits(JAVA["raw_membership_bits"], (np.prod(shape) * 4,)).view(np.int32),
    )
    assert raw.objective.view(np.uint32) == np.asarray(JAVA["raw_objective_bits"], dtype=np.uint32)
    assert not np.array_equal(raw.memberships, result.memberships)


def test_top_four_ties_zero_energy_and_ineligible_prior_match_source_order() -> None:
    result = _run()
    assert result.best_class_indices[9, 9, 9].tolist() == JAVA["tie_best"]
    assert result.best_class_indices[8, 8, 8].tolist() == JAVA["zero_energy_best"]
    assert result.best_class_indices[7, 7, 7].tolist() == JAVA["near_zero_best"]
    assert result.best_class_indices[10, 10, 10].tolist() == JAVA["no_eligible_best"]
    assert result.memberships[8, 8, 8, 0] == np.float32(1.0)
    assert result.memberships[7, 7, 7, 0] == np.float32(1.0)
    assert result.denominator[10, 10, 10] == np.float32(0.0)
    assert np.array_equal(result.memberships[10, 10, 10], np.zeros(4, np.float32))


def test_six_face_candidates_affect_smoothing_and_borders_remain_untouched() -> None:
    baseline = _run()
    original_memberships = _nifti("memberships")
    original_best = _nifti("best").astype(np.uint8)
    target = (1, 1, 1)
    target_class = int(baseline.best_class_indices[target][0])
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    changes = []
    for dx, dy, dz in offsets:
        changed = original_memberships.copy()
        neighbor = (target[0] + dx, target[1] + dy, target[2] + dz)
        slot = int(np.flatnonzero(original_best[neighbor] != target_class)[0])
        changed[neighbor + (slot,)] = np.float32(changed[neighbor + (slot,)] + 0.17)
        altered = _run(overrides={"initial_memberships": changed})
        changes.append(altered.denominator[target] != baseline.denominator[target])
    assert all(changes)

    for axis in range(3):
        for index in (0, -1):
            slices = [slice(None)] * 4
            slices[axis] = index
            np.testing.assert_array_equal(baseline.memberships[tuple(slices)], original_memberships[tuple(slices)])
            np.testing.assert_array_equal(baseline.best_class_indices[tuple(slices)], original_best[tuple(slices)])
