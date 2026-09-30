from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import compute_lesion_toads_auxiliary_centroids


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_special_centroids"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "LesionAuxiliaryCentroidsReferenceHarness.java"
HARNESS_SHA256 = "0dca4a600bf029ec9695cbbfbf743ea35d8055b324f2b359e8e03d3f500c4594"
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
INITIAL = np.asarray(
    [[0.15, 0.55, 1.05, 1.65, 2.0], [0.30, 0.70, 1.15, 1.80, 2.1]], dtype=np.float32
)
PRIOR = np.asarray(
    [[0.20, 0.60, 0.95, 1.40, 2.2], [0.35, 0.75, 1.10, 1.55, 2.3]], dtype=np.float32
)
LESION_CENTROID = np.asarray([0.85, 0.95], dtype=np.float32)
BLACK_CENTROID = np.asarray([0.45, 0.52], dtype=np.float32)
LESION_EXPONENTS = np.asarray([1.2, 0.8], dtype=np.float32)
LINKS = np.asarray([[0, 1, 1, 2, 3], [0, 0, 1, 1, 2]], dtype=np.int16)
TYPES = ("object", "lesion", "blackhole", "object", "object")


def _bits(values: np.ndarray) -> list[list[int]] | list[int]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def _inputs() -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
]:
    image = nib.load(str(ROOT / "special_centroid_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 4, 3, 12)
    return (
        frames[..., 0].astype(np.int16),
        frames[..., 1:6],
        frames[..., 6:8],
        frames[..., 8:10],
        frames[..., 10],
        frames[..., 11],
        frames,
    )


def test_lesion_and_blackhole_centroid_routing_match_java_powertable() -> None:
    assert REFERENCE["source_sha256"] == "77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0"
    assert REFERENCE["powertable_sha256"] == POWER_TABLE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert all(value > 0 for value in REFERENCE["route_counts"])
    assert all(value > 0 for value in REFERENCE["tie_counts"])
    labels, memberships, images, fields, lesion, black, frames = _inputs()
    assert set(np.unique(labels)) == {0, 1, 2, 3}  # class 4 has no hard-label support
    # Assert fixture locations survive serialization and drive the documented
    # source tie precedence and a zero-distance lesion-winner branch.
    assert labels[0, 0, 1] == 1
    assert memberships[0, 0, 1, 1] == lesion[0, 0, 1] > black[0, 0, 1]
    assert labels[0, 0, 2] == 2
    assert lesion[0, 0, 2] == black[0, 0, 2] > memberships[0, 0, 2, 2]
    assert labels[1, 0, 1] == 2
    assert memberships[1, 0, 1, 2] == lesion[1, 0, 1] == black[1, 0, 1]
    assert labels[0, 1, 0] == 2 and black[0, 1, 0] > lesion[0, 1, 0]
    assert labels[1, 0, 0] == 1 and lesion[1, 0, 0] > memberships[1, 0, 0, 1]
    assert fields[1, 0, 0, 0] == np.float32(1.0)
    assert images[1, 0, 0, 0] == LESION_CENTROID[0]

    for mode in ("none", "linked", "damped", "prior", "fallback"):
        for field_key, field_values in (("no_field", None), ("with_field", fields)):
            kwargs: dict[str, object] = {"fields": field_values, "mode": mode, "smoothness": 0.35}
            if mode == "linked":
                kwargs["centroid_links"] = LINKS
            if mode == "prior":
                kwargs["prior_centroids"] = PRIOR
            actual = compute_lesion_toads_auxiliary_centroids(
                images, memberships, labels, INITIAL, lesion, black,
                LESION_CENTROID, BLACK_CENTROID, TYPES, LESION_EXPONENTS, **kwargs,
            )
            expected = REFERENCE[f"{mode}_{field_key}"]
            assert _bits(actual.class_centroids) == expected["class"]
            assert _bits(actual.lesion_centroids) == expected["lesion"]
            assert _bits(actual.black_hole_centroids) == expected["blackhole"]


def test_special_centroids_reject_malformed_maps_and_parameters() -> None:
    labels, memberships, images, fields, lesion, black, _frames = _inputs()
    args = (images, memberships, labels, INITIAL, lesion, black,
            LESION_CENTROID, BLACK_CENTROID, TYPES, LESION_EXPONENTS)
    with pytest.raises(ValueError, match="match the image grid"):
        compute_lesion_toads_auxiliary_centroids(
            *args[:5], lesion[..., :2], *args[6:]
        )
    with pytest.raises(ValueError, match="lesion_power_exponents"):
        compute_lesion_toads_auxiliary_centroids(
            *args[:-1], np.asarray([1.2], dtype=np.float32)
        )
    with pytest.raises(ValueError, match="fields"):
        compute_lesion_toads_auxiliary_centroids(
            *args, fields=fields[..., :1]
        )
    with pytest.raises(ValueError, match="unsupported"):
        compute_lesion_toads_auxiliary_centroids(
            images, memberships, labels, INITIAL, lesion, black,
            LESION_CENTROID, BLACK_CENTROID,
            ("object", "optimized_gm", "blackhole", "object", "object"),
            LESION_EXPONENTS,
            optimized_power_exponents=LESION_EXPONENTS,
        )
