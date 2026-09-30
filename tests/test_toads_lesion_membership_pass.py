"""Selected-Java comparison for the fixed-input base-TOADS LESION branch."""

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.toads_lesion_membership_pass import toads_lesion_membership_pass


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_lesion_membership_pass"
SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
MODALITY_WEIGHTS = np.asarray([[0.75, 0.35], [1.125, 0.625]], dtype=np.float32)
IMAGE_VARIANCES = np.asarray([0.8, 1.25], dtype=np.float32)
VARIANCE_WEIGHTS = np.asarray([[0.8, 1.1], [1.2, 0.7]], dtype=np.float32)
CENTROIDS = np.asarray([[0.25, 0.125], [-0.125, -0.25]], dtype=np.float32)
LESION_CENTROID = np.asarray([0.125, -0.25], dtype=np.float32)
SMOOTH_PRIORS = np.asarray([0.6, 0.9], dtype=np.float32)
ATLAS_PRIORS = np.asarray([0.21, 0.34], dtype=np.float32)
SIMILARITY = np.asarray([[0.0, 0.4], [0.35, 0.0]], dtype=np.float32)
CASE_SETTINGS = {
    "object_lesion_field_bias_smooth_prior": (False, True, True, 0.18, 0.31, False),
    "lesion_object_direct_no_field_no_bias": (True, False, False, 0.0, 0.0, True),
    "object_lesion_field_smoothing_only": (False, True, False, 0.14, 0.0, True),
    "lesion_object_bias_prior_only": (True, False, True, 0.0, 0.22, False),
}


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _expected(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def _case_inputs(inputs: np.ndarray, reverse: bool, all_direct: bool):
    order = np.asarray([1, 0] if reverse else [0, 1])
    types = ("LESION", "OBJECT") if reverse else ("OBJECT", "LESION")
    initial = inputs[..., 6:8][..., order]
    shapes = inputs[..., 9:11][..., order]
    lesion_index = types.index("LESION")
    object_index = types.index("OBJECT")
    segmentation = np.where(inputs[..., 11] == np.float32(1.0), lesion_index, object_index).astype(np.int8)
    closest = inputs[..., 12:16].astype(np.int8)
    if reverse:
        closest = np.where(closest >= 0, 1 - closest, closest).astype(np.int8)
    if all_direct:
        closest[...] = np.asarray([0, 1, -1, -1], dtype=np.int8)
    return types, initial, inputs[..., 8], segmentation, inputs[..., :2], shapes, closest


def test_object_lesion_pass_matches_source_excerpt_bitwise() -> None:
    image = nib.load(str(ROOT / "inputs_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    inputs = np.asarray(image.dataobj, dtype=np.float32)
    assert inputs.shape == (5, 4, 5, 16)
    reference = json.loads((ROOT / "java_reference.json").read_text())
    digest = hashlib.sha256(inputs.astype("<f4", copy=False).ravel(order="C").view("<u4").tobytes()).hexdigest()
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["input_sha256"] == digest
    assert reference["shape"] == [5, 4, 5]
    assert reference["channels"] == 16

    for case in reference["cases"]:
        reverse, use_field, use_bias, smoothing, prior_coefficient, all_direct = CASE_SETTINGS[case["name"]]
        types, initial, lesion_initial, segmentation, images, shapes, closest = _case_inputs(
            inputs, reverse, all_direct
        )
        order = np.asarray([1, 0] if reverse else [0, 1])
        result = toads_lesion_membership_pass(
            images,
            initial,
            lesion_initial,
            segmentation,
            types,
            CENTROIDS[:, order],
            LESION_CENTROID,
            MODALITY_WEIGHTS,
            VARIANCE_WEIGHTS[:, order],
            IMAGE_VARIANCES,
            shapes,
            closest,
            SMOOTH_PRIORS[order],
            SIMILARITY[np.ix_(order, order)],
            ATLAS_PRIORS[order],
            smoothing,
            0.57,
            prior_coefficient,
            intensity_field=inputs[..., 2:4] if use_field else None,
            psf_bias=inputs[..., 4:6] if use_bias else None,
        )
        for name, actual in (
            ("class_data_terms_bits", result.class_data_terms),
            ("lesion_data_terms_bits", result.lesion_data_terms),
            ("class_energy_bits", result.class_energy),
            ("memberships_bits", result.memberships),
            ("lesion_membership_bits", result.lesion_membership),
            ("denominator_bits", result.denominator),
        ):
            expected = _expected(case[name], actual.shape)
            np.testing.assert_array_equal(actual.view(np.int32), expected.view(np.int32), err_msg=f"{case['name']}:{name}")
        assert _bits(np.asarray([result.objective])) == [case["objective_bits"]]
        assert _bits(result.class_data_terms) == case["class_data_terms_bits"]

        # Each source pass leaves all six boundary faces untouched.
        for axis in range(3):
            for edge in (0, -1):
                key = [slice(None)] * 3
                key[axis] = edge
                key = tuple(key)
                np.testing.assert_array_equal(result.memberships[key], initial[key])
                np.testing.assert_array_equal(result.lesion_membership[key], lesion_initial[key])
                assert not np.any(result.denominator[key])

    zero_case = next(case for case in reference["cases"] if case["name"] == "lesion_object_direct_no_field_no_bias")
    zero_terms = _expected(zero_case["lesion_data_terms_bits"], (5, 4, 5, 2))
    class_terms = _expected(zero_case["class_data_terms_bits"], (5, 4, 5, 2))
    assert zero_terms[2, 1, 2, 0] == 0.0  # image equals lesion centroid on assigned voxel
    assert class_terms[3, 1, 2, 1] == 0.0  # image equals OBJECT centroid


def test_lesion_pass_rejects_inactive_prior_only_shape() -> None:
    inputs = np.asarray(nib.load(str(ROOT / "inputs_ras.nii.gz")).dataobj, dtype=np.float32)
    types, initial, lesion, segmentation, images, shapes, closest = _case_inputs(inputs, False, False)
    # Force one class to PRIOR at this interior voxel and suppress its shape.
    closest[2, 1, 2] = (-1, -1, -1, -1)
    shapes[2, 1, 2, 0] = 0.0
    with np.testing.assert_raises_regex(ValueError, "all interior classes must have active shapes"):
        toads_lesion_membership_pass(
            images, initial, lesion, segmentation, types, CENTROIDS, LESION_CENTROID,
            MODALITY_WEIGHTS, VARIANCE_WEIGHTS, IMAGE_VARIANCES, shapes, closest,
            SMOOTH_PRIORS, SIMILARITY, ATLAS_PRIORS, 0.1, 0.57, 0.2,
        )
