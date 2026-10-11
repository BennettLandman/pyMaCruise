from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_standard13_voxel import (
    STANDARD13_NAMES,
    STANDARD13_TYPES,
    STANDARD13_CLASS_SMOOTHING,
    STANDARD13_PRIOR_COEFFICIENT,
    STANDARD13_PRIOR_SCALE,
    compute_lesion_toads_standard13_voxel,
)

ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_standard13_voxel"
WRAPPER_EXCERPT = ROOT / "wrapper_initial_membership_source_excerpt.java.txt"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
WRAPPER_EXCERPT_SHA = "ffe1304226cad722ea12c12654a31889c2f988870e64105e7f2a84aecb43fd6c"


def _state() -> dict[str, np.ndarray | float | int]:
    image = np.asarray(nib.load(str(ROOT / "inputs_ras.nii.gz")).dataobj, dtype=np.float32)[1, 1, 1]
    k = np.arange(13, dtype=np.float32)
    centers = np.stack((np.float32(0.2) + np.float32(0.04) * k,
                        np.float32(0.3) + np.float32(0.03) * k))
    deviations = np.stack((np.float32(0.3) + np.float32(0.01) * k,
                           np.float32(0.4) + np.float32(0.005) * k))
    return {
        "image_values": image,
        "class_centroids": centers,
        "class_stddev": deviations,
        "lesion_centroid": np.asarray([0.5, 0.65], dtype=np.float32),
        "lesion_stddev": np.asarray([0.25, 0.3], dtype=np.float32),
        "channel_weights": np.asarray([1.0, 0.75], dtype=np.float32),
        "shape_values": np.float32(0.4) + np.float32(0.02) * k,
        "shape_priors": np.float32(0.5) + np.float32(0.05) * k,
        "smoothing_priors": np.float32(0.8) + np.float32(0.025) * k,
        "direct_relations": np.asarray([True, False, True, False, True, False, True,
                                         False, True, False, True, False, True]),
        "class_memberships": np.float32(0.01) + np.float32(0.001) * k,
        "lesion_membership": np.float32(0.125),
        "black_hole_membership": np.float32(0.0),
        "outlier_membership": np.float32(0.0),
        "lesion_relation": np.float32(0.0),
        "six_face_class_memberships": np.asarray(
            [[np.float32(0.02) + np.float32(0.001) * (face * 13 + index)
              for index in range(13)] for face in range(6)], dtype=np.float32),
        "six_face_lesion_memberships": np.asarray(
            [0.05, 0.1, 0.15, 0.2, 0.25, 0.3], dtype=np.float32),
        "six_face_black_hole_memberships": np.asarray(
            [0.0, 0.01, 0.02, 0.03, 0.04, 0.05], dtype=np.float32),
        "outlier_energy": np.float32(0.1),
        "lesion_power_exponent": np.float32(1.2),
        "ventricle_distance": 1,
        "gray_matter_distance": -1,
        "brainstem_distance": 10,
    }


def _bits(value: np.ndarray | np.float32) -> list[int] | int:
    values = np.asarray(value, dtype=np.float32)
    if values.ndim == 0:
        return int(values.view(np.uint32))
    return values.view(np.uint32).reshape(-1).astype(np.uint64).tolist()


def test_standard13_profile_composes_selected_initial_membership_voxel() -> None:
    image = nib.load(str(ROOT / "inputs_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    image_values = np.asarray(image.dataobj, dtype=np.float32)
    image_sha = hashlib.sha256(image_values.astype("<f4", copy=False).tobytes()).hexdigest()
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["wrapper_source_sha256"] == WRAPPER_SHA
    wrapper_excerpt = WRAPPER_EXCERPT.read_bytes()
    assert hashlib.sha256(wrapper_excerpt).hexdigest() == WRAPPER_EXCERPT_SHA
    wrapper_text = wrapper_excerpt.decode("utf-8")
    wrapper_calls = [
        "segmentation.computeRelations();",
        "segmentation.computeVariances();",
        "segmentation.computeChannelWeights();",
        "energy = segmentation.computeMemberships();",
        "segmentation.setGMflag();",
    ]
    assert [wrapper_text.index(call) for call in wrapper_calls] == sorted(
        wrapper_text.index(call) for call in wrapper_calls
    )
    assert reference["image_sha256"] == image_sha
    assert reference["class_names"] == list(STANDARD13_NAMES)
    assert reference["class_types"] == list(STANDARD13_TYPES)
    assert reference["profile_bits"] == {
        "prior_coefficient": int(np.asarray(STANDARD13_PRIOR_COEFFICIENT).view(np.uint32)),
        "prior_scale": int(np.asarray(STANDARD13_PRIOR_SCALE).view(np.uint32)),
        "class_smoothing": int(np.asarray(STANDARD13_CLASS_SMOOTHING).view(np.uint32)),
    }

    state = _state()
    result = compute_lesion_toads_standard13_voxel(**state)
    actual = {
        "class_membership_bits": _bits(result.class_memberships),
        "lesion_membership_bits": _bits(result.lesion_membership),
        "black_hole_membership_bits": _bits(result.black_hole_membership),
        "outlier_membership_bits": _bits(result.outlier_membership),
        "lesion_relation_bits": _bits(result.lesion_relation),
        "class_energy_bits": _bits(result.class_energies),
        "lesion_energy_bits": _bits(result.lesion_energy),
        "black_hole_energy_bits": _bits(result.black_hole_energy),
        "denominator_bits": _bits(result.denominator),
        "objective_bits": _bits(result.objective),
        "class_similarity_bits": _bits(result.class_similarity),
        "lesion_similarity_bits": _bits(result.lesion_similarity),
    }
    for name, values in actual.items():
        assert values == reference[name], name


def test_standard13_zero_shape_clears_skipped_class_membership() -> None:
    state = _state()
    shapes = np.asarray(state["shape_values"], dtype=np.float32).copy()
    shapes[8] = 0.0
    direct = np.asarray(state["direct_relations"], dtype=np.bool_).copy()
    direct[8] = False
    state["shape_values"] = shapes
    state["direct_relations"] = direct
    result = compute_lesion_toads_standard13_voxel(**state)
    assert result.class_memberships[8] == np.float32(0.0)
    assert result.class_energies[8] == 0.0


@pytest.mark.parametrize("direct_relation", [False, True])
def test_zero_shape_leaves_lesion_auxiliary_state_like_java_guard(direct_relation: bool) -> None:
    state = _state()
    shapes = np.asarray(state["shape_values"], dtype=np.float32).copy()
    shapes[10] = np.float32(0.0)
    direct = np.asarray(state["direct_relations"], dtype=np.bool_).copy()
    direct[10] = direct_relation
    state["shape_values"] = shapes
    state["direct_relations"] = direct
    state["lesion_membership"] = np.float32(0.25)
    state["black_hole_membership"] = np.float32(0.5)
    state["lesion_relation"] = np.float32(0.75)

    result = compute_lesion_toads_standard13_voxel(**state)

    assert result.class_memberships[10] == np.float32(0.0)
    assert result.class_energies[10] == np.float32(0.0)
    assert result.lesion_membership == np.float32(np.float32(0.25) / result.denominator)
    assert result.black_hole_membership == np.float32(0.5)
    assert result.lesion_relation == np.float32(0.75)


def test_source_zero_variance_mode_is_valid_when_current_statistics_are_nonzero() -> None:
    state = _state()
    source = compute_lesion_toads_standard13_voxel(
        **state,
        gm_flag=True,
        source_zero_variance_nan=True,
    )
    finite = compute_lesion_toads_standard13_voxel(
        **state,
        gm_flag=True,
        source_zero_variance_nan=False,
    )
    for name in (
        "class_memberships",
        "lesion_membership",
        "black_hole_membership",
        "outlier_membership",
        "lesion_relation",
        "class_energies",
        "lesion_energy",
        "black_hole_energy",
        "denominator",
    ):
        np.testing.assert_array_equal(getattr(source, name), getattr(finite, name))
    assert source.objective.view(np.uint32) == finite.objective.view(np.uint32)


def test_standard13_voxel_field_application_matches_unit_and_changes_residuals() -> None:
    state = _state()
    uncorrected = compute_lesion_toads_standard13_voxel(**state)
    identity = compute_lesion_toads_standard13_voxel(
        **state, intensity_field=np.ones(2, dtype=np.float32)
    )
    np.testing.assert_array_equal(identity.class_memberships, uncorrected.class_memberships)
    np.testing.assert_array_equal(identity.class_energies, uncorrected.class_energies)
    assert identity.objective.view(np.uint32) == uncorrected.objective.view(np.uint32)

    corrected = compute_lesion_toads_standard13_voxel(
        **state, intensity_field=np.asarray([1.1, 0.9], dtype=np.float32)
    )
    assert corrected.class_energies[3].view(np.uint32) != uncorrected.class_energies[3].view(np.uint32)
    assert corrected.class_memberships[3].view(np.uint32) != uncorrected.class_memberships[3].view(np.uint32)


def test_standard13_voxel_rejects_invalid_correction_field() -> None:
    state = _state()
    for field in (
        np.ones(3, dtype=np.float32),
        np.asarray([1.0, np.nan], dtype=np.float32),
    ):
        try:
            compute_lesion_toads_standard13_voxel(**state, intensity_field=field)
        except ValueError as error:
            assert "intensity_field" in str(error)
        else:  # pragma: no cover - assertion branch
            raise AssertionError("invalid intensity field was accepted")
