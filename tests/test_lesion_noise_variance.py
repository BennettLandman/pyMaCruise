from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import compute_lesion_toads_noise_standard_deviation


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_variances"
REFERENCE = json.loads((ROOT / "noise_java_reference.json").read_text())
CENTROIDS = np.asarray(
    [[0.15, 0.55, 1.05, 1.65], [0.30, 0.70, 1.15, 1.80]], dtype=np.float32
)


def _bits(values: np.ndarray) -> list[list[int]]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def test_pooled_noise_standard_deviation_matches_java() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    image = nib.load(str(ROOT / "variance_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    labels = frames[..., 0].astype(np.int16)
    memberships = frames[..., 1:5]
    values = frames[..., 5:7]
    fields = frames[..., 7:9]
    no_field = compute_lesion_toads_noise_standard_deviation(
        values, memberships, labels, CENTROIDS
    )
    with_field = compute_lesion_toads_noise_standard_deviation(
        values, memberships, labels, CENTROIDS, fields=fields
    )
    assert _bits(no_field) == REFERENCE["no_field"]
    assert _bits(with_field) == REFERENCE["with_field"]
    for result in (no_field, with_field):
        np.testing.assert_array_equal(result, np.broadcast_to(result[:, :1], result.shape))
    altered = values.copy()
    altered[labels == 0] = np.float32(100000.0)
    np.testing.assert_array_equal(
        compute_lesion_toads_noise_standard_deviation(
            altered, memberships, labels, CENTROIDS
        ),
        no_field,
    )

    empty_support = compute_lesion_toads_noise_standard_deviation(
        values, memberships, np.zeros_like(labels), CENTROIDS
    )
    assert _bits(empty_support) == REFERENCE["empty_support"]
    np.testing.assert_array_equal(empty_support, np.zeros_like(empty_support))


def test_noise_variance_rejects_bad_centroid_shape() -> None:
    volume = np.zeros((3, 3, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="class_centroids"):
        compute_lesion_toads_noise_standard_deviation(
            volume[..., None], np.ones((*volume.shape, 2)),
            np.zeros(volume.shape, dtype=np.int16), np.zeros((2, 2)),
        )
