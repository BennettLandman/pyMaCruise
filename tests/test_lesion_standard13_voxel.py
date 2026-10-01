from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
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
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


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


def test_standard13_profile_rejects_nonpositive_shape_for_active_class_loop() -> None:
    state = _state()
    shapes = np.asarray(state["shape_values"], dtype=np.float32).copy()
    shapes[8] = 0.0
    state["shape_values"] = shapes
    with np.testing.assert_raises_regex(ValueError, "all 13 classes must have positive"):
        compute_lesion_toads_standard13_voxel(**state)
