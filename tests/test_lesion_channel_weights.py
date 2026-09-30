from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import compute_lesion_toads_channel_weights


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_channel_weights"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
TYPES = ("object", "outlier", "lesion", "blackhole", "mask", "outmask")
CENTROIDS = np.asarray(
    [[0.05, 0.10, 0.15, 0.20, 0.0, 0.0], [0.10, 0.15, 0.20, 0.25, 0.0, 0.0]],
    dtype=np.float32,
)
LESION_CENTROID = np.asarray([0.03, 0.04], dtype=np.float32)
BLACK_HOLE_CENTROID = np.asarray([0.02, 0.02], dtype=np.float32)


def _bits(values: np.ndarray) -> list[int]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def test_diagonal_channel_weights_match_selected_java() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    image = nib.load(str(ROOT / "channel_weight_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 4, 3, 13)
    memberships = frames[..., :6]
    values = frames[..., 6:8]
    fields = frames[..., 8:10]
    lesion = frames[..., 10]
    black_hole = frames[..., 11]
    outlier = frames[..., 12]
    inputs_before = tuple(array.copy() for array in (memberships, values, fields, lesion, black_hole, outlier))

    no_field = compute_lesion_toads_channel_weights(
        values, memberships, CENTROIDS, lesion, black_hole, outlier,
        LESION_CENTROID, BLACK_HOLE_CENTROID, TYPES, outlier_distance=2.75,
    )
    with_field = compute_lesion_toads_channel_weights(
        values, memberships, CENTROIDS, lesion, black_hole, outlier,
        LESION_CENTROID, BLACK_HOLE_CENTROID, TYPES, outlier_distance=2.75,
        fields=fields,
    )
    uniform = compute_lesion_toads_channel_weights(
        values, memberships, CENTROIDS, lesion, black_hole, outlier,
        LESION_CENTROID, BLACK_HOLE_CENTROID, TYPES, outlier_distance=2.75,
        norm_type="Diagonal",
    )
    zero_energy = compute_lesion_toads_channel_weights(
        values, memberships[..., :2], CENTROIDS[:, :2], lesion, black_hole,
        outlier, LESION_CENTROID, BLACK_HOLE_CENTROID, ("mask", "outmask"),
        outlier_distance=2.75,
    )
    assert _bits(no_field) == REFERENCE["diagonal_no_field"]
    assert _bits(with_field) == REFERENCE["diagonal_with_field"]
    assert _bits(uniform) == REFERENCE["uniform"]
    assert _bits(zero_energy) == REFERENCE["zero_energy"]
    np.testing.assert_array_equal(uniform, np.asarray([0.5, 0.5], dtype=np.float32))
    np.testing.assert_array_equal(zero_energy, uniform)
    for actual, before in zip((memberships, values, fields, lesion, black_hole, outlier), inputs_before):
        np.testing.assert_array_equal(actual, before)


def test_channel_weights_reject_malformed_maps() -> None:
    image = np.zeros((3, 3, 3, 1), dtype=np.float32)
    memberships = np.ones((3, 3, 3, 2), dtype=np.float32)
    maps = np.zeros((3, 3, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="object_types"):
        compute_lesion_toads_channel_weights(
            image, memberships, np.zeros((1, 2)), maps, maps, maps,
            np.zeros(1), np.zeros(1), ("object",), outlier_distance=2.0,
        )
