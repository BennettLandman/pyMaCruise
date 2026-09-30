from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import compute_lesion_toads_centroids


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_mixed_centroids"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "LesionMixedCentroidsReferenceHarness.java"
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
HARNESS_SHA256 = "d80f4b7b93ba547e718ea8fa5dc37e8c5a5b916a5bc20b7d7ebd28059f9c7dcf"
TYPES = (
    "object", "mask", "outmask", "outlier", "optimized_csf", "optimized_gm",
    "optimized_wm", "lesion", "blackhole", "optimized_lesion",
    "optimized_blackhole", "rayleigh", "object",
)
INITIAL = np.asarray(
    [[np.float32(np.float32(0.15) + np.float32(np.float32(0.17) * k) + np.float32(np.float32(0.15) * c)) for k in range(13)] for c in range(2)],
    dtype=np.float32,
)
PRIOR = np.asarray(
    [[np.float32(np.float32(0.2) + np.float32(np.float32(0.13) * k) + np.float32(np.float32(0.1) * c)) for k in range(13)] for c in range(2)],
    dtype=np.float32,
)
TISSUE_POWER = np.asarray(
    [[np.float32(np.float32(0.75) + np.float32(np.float32(0.11) * ((k + 2 * c) % 7))) for k in range(13)] for c in range(2)],
    dtype=np.float32,
)
LINKS = np.asarray([[(k + 2 * c) % 5 for k in range(13)] for c in range(2)], dtype=np.int16)
LESION_CENTROIDS = np.asarray([0.85, 0.95], dtype=np.float32)
BLACK_CENTROIDS = np.asarray([0.45, 0.52], dtype=np.float32)


def _bits(values: np.ndarray) -> list[list[int]] | list[int]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def _inputs() -> tuple[np.ndarray, ...]:
    nii = nib.load(str(ROOT / "mixed_centroid_inputs.nii.gz"))
    assert nib.aff2axcodes(nii.affine) == ("R", "A", "S")
    frames = np.asarray(nii.dataobj, dtype=np.float32)
    assert frames.shape == (6, 5, 4, 21)
    return (
        frames[..., 0].astype(np.int16), frames[..., 1:14], frames[..., 14:16],
        frames[..., 16:18], frames[..., 18], frames[..., 19], frames[..., 20],
    )


def test_mixed_centroid_modes_match_selected_java_bits() -> None:
    assert REFERENCE["source_sha256"] == "77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0"
    assert REFERENCE["powertable_sha256"] == POWER_TABLE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    labels, memberships, images, fields, outlier, lesion, black = _inputs()
    assert set(np.unique(labels)) == set(range(12))  # class 12 is intentionally absent
    assert len(set(TYPES)) == 12
    for mode in ("none", "linked", "damped", "prior", "fallback"):
        for key, field_values in (("no_field", None), ("with_field", fields)):
            kwargs: dict[str, object] = {
                "fields": field_values, "mode": mode, "smoothness": 0.35,
                "outlier_membership": outlier, "lesion_membership": lesion,
                "black_hole_membership": black, "lesion_centroids": LESION_CENTROIDS,
                "black_hole_centroids": BLACK_CENTROIDS,
                "tissue_power_exponents": TISSUE_POWER,
                "modalities": ("T1_MPRAGE", "FLAIR"),
                "lesion_power_exponents": np.asarray([1.2, 0.8], dtype=np.float32),
                "optimized_auxiliary_power_exponents": np.asarray([1.25, 0.75], dtype=np.float32),
            }
            if mode == "linked":
                kwargs["centroid_links"] = LINKS
            if mode == "prior":
                kwargs["prior_centroids"] = PRIOR
            actual = compute_lesion_toads_centroids(
                images, memberships, labels, INITIAL, TYPES, **kwargs
            )
            expected = REFERENCE[f"{mode}_{key}"]
            assert _bits(actual.class_centroids) == expected["class"]
            assert actual.lesion_centroids is not None
            assert actual.black_hole_centroids is not None
            assert _bits(actual.lesion_centroids) == expected["lesion"]
            assert _bits(actual.black_hole_centroids) == expected["blackhole"]


def test_mixed_centroids_validate_category_specific_inputs() -> None:
    labels, memberships, images, _fields, outlier, lesion, black = _inputs()
    kwargs: dict[str, object] = {
        "outlier_membership": outlier, "lesion_membership": lesion,
        "black_hole_membership": black, "lesion_centroids": LESION_CENTROIDS,
        "black_hole_centroids": BLACK_CENTROIDS,
        "tissue_power_exponents": TISSUE_POWER, "modalities": ("T1_MPRAGE", "FLAIR"),
        "lesion_power_exponents": np.asarray([1.2, 0.8], dtype=np.float32),
        "optimized_auxiliary_power_exponents": np.asarray([1.25, 0.75], dtype=np.float32),
    }
    with pytest.raises(ValueError, match="outlier map"):
        compute_lesion_toads_centroids(
            images, memberships, labels, INITIAL, TYPES,
            outlier_membership=None, **{k: v for k, v in kwargs.items() if k != "outlier_membership"},
        )
    with pytest.raises(ValueError, match="auxiliary class types"):
        compute_lesion_toads_centroids(
            images, memberships, labels, INITIAL, TYPES,
            **{k: v for k, v in kwargs.items() if k not in {"lesion_membership", "black_hole_membership"}},
        )
    with pytest.raises(ValueError, match="tissue_power_exponents"):
        compute_lesion_toads_centroids(
            images, memberships, labels, INITIAL, TYPES,
            **{k: v for k, v in kwargs.items() if k != "tissue_power_exponents"},
        )
    with pytest.raises(ValueError, match="optimized_auxiliary_power_exponents"):
        compute_lesion_toads_centroids(
            images, memberships, labels, INITIAL, TYPES,
            **{k: v for k, v in kwargs.items() if k != "optimized_auxiliary_power_exponents"},
        )
