from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_statistics_refresh import refresh_lesion_toads_statistics


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_statistics_refresh"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "LesionStatisticsRefreshReferenceHarness.java"
TYPES = ("object", "optimized_gm", "lesion", "blackhole", "object")
INITIAL = np.asarray(
    [[
        np.float32(
            np.float32(np.float32(0.2) + np.float32(np.float32(0.15) * k))
            + np.float32(np.float32(0.11) * c)
        )
        for k in range(5)
    ] for c in range(2)],
    dtype=np.float32,
)
INITIAL_LESION = np.asarray([0.85, 0.95], dtype=np.float32)
INITIAL_BLACK = np.asarray([0.45, 0.52], dtype=np.float32)
TISSUE_POWER = np.asarray([[0.9, 0.9, 0.9, 0.9, 0.9], [1.05] * 5], dtype=np.float32)
LESION_POWER = np.asarray([1.2, 0.8], dtype=np.float32)


def _bits(values: np.ndarray) -> list[int] | list[list[int]]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def _inputs() -> tuple[np.ndarray, ...]:
    nii = nib.load(str(ROOT / "statistics_refresh_inputs.nii.gz"))
    assert nib.aff2axcodes(nii.affine) == ("R", "A", "S")
    frames = np.asarray(nii.dataobj, dtype=np.float32)
    assert frames.shape == (6, 5, 4, 13)
    return (
        frames[..., 1:6], frames[..., 0].astype(np.int16), frames[..., 6:8],
        frames[..., 8:10], frames[..., 10], frames[..., 11], frames[..., 12],
    )


def test_statistics_refresh_matches_selected_java_wrapper_order() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert REFERENCE["wrapper_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    power_table = Path(__file__).parent / "java_reference" / "upstream" / "edu" / "jhmi" / "rad" / "medic" / "utilities" / "PowerTable.java"
    assert hashlib.sha256(power_table.read_bytes()).hexdigest() == REFERENCE["powertable_sha256"]
    assert REFERENCE["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    memberships, labels, images, fields, outlier, lesion, black = _inputs()
    assert set(np.unique(labels)) == {0, 1, 2, 3}  # class 4 has no hard labels
    assert memberships.shape[-1] == len(TYPES)
    saved = tuple(a.copy() for a in (memberships, labels, images, fields, outlier, lesion, black))

    centroid_options = {
        "mode": "none",
        "tissue_power_exponents": TISSUE_POWER,
        "modalities": ("T1_MPRAGE", "FLAIR"),
        "lesion_power_exponents": LESION_POWER,
    }
    for key, field_values in (("no_field", None), ("with_field", fields)):
        actual = refresh_lesion_toads_statistics(
            images,
            memberships,
            labels,
            INITIAL,
            lesion,
            black,
            outlier,
            INITIAL_LESION,
            INITIAL_BLACK,
            TYPES,
            fields=field_values,
            centroid_options=centroid_options,
            channel_weight_options={"outlier_distance": 2.75},
        )
        expected = REFERENCE[key]
        assert _bits(actual.centroids.class_centroids) == expected["class"]
        assert actual.centroids.lesion_centroids is not None
        assert actual.centroids.black_hole_centroids is not None
        assert _bits(np.stack((actual.centroids.lesion_centroids, actual.centroids.black_hole_centroids))) == expected["aux_centroids"]
        assert _bits(actual.variances.class_standard_deviation.reshape(-1)) == expected["class_std"]
        assert _bits(actual.variances.lesion_standard_deviation) == expected["lesion_std"]
        assert _bits(actual.variances.black_hole_standard_deviation) == expected["blackhole_std"]
        assert _bits(actual.channel_weights) == expected["weights"]
        assert actual.centroids.class_centroids[:, 4].tolist() == [0.0, 0.0]
        assert actual.variances.class_standard_deviation[:, 4].tolist() == [0.0, 0.0]
        assert np.unique(actual.channel_weights).size == 2

    fallback = refresh_lesion_toads_statistics(
        images, memberships, labels, INITIAL, lesion, black, outlier,
        INITIAL_LESION, INITIAL_BLACK, TYPES,
        centroid_options=centroid_options,
        channel_weight_options={"outlier_distance": 2.75, "norm_type": "Diagonal"},
    )
    assert _bits(fallback.channel_weights) == REFERENCE["uniform_weight_fallback"]

    for actual, before in zip((memberships, labels, images, fields, outlier, lesion, black), saved):
        np.testing.assert_array_equal(actual, before)
