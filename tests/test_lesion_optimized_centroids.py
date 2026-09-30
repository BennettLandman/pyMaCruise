from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import compute_lesion_toads_basic_centroids


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_optimized_centroids"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "LesionOptimizedCentroidsReferenceHarness.java"
HARNESS_SHA256 = "414d9b40320305523c14dc3114c66ae687a271f104f2dce6cc3ef31c9e346dbe"
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
INITIAL = np.asarray(
    [[0.15, 0.55, 1.05, 1.65, 2.0], [0.30, 0.70, 1.15, 1.80, 2.1]], dtype=np.float32
)
PRIOR = np.asarray(
    [[0.20, 0.60, 0.95, 1.40, 2.2], [0.35, 0.75, 1.10, 1.55, 2.3]], dtype=np.float32
)
EXPONENTS = np.asarray(
    [[0.0, 1.35, 0.8, 1.2, 0.0], [0.0, 1.1, 1.25, 0.75, 0.0]], dtype=np.float32
)
LINKS = np.asarray([[0, 1, 1, 2, 3], [0, 0, 1, 1, 2]], dtype=np.int16)
TYPES = ("object", "optimized_csf", "optimized_gm", "optimized_wm", "object")
MODALITIES = ("T1_MPRAGE", "FLAIR")


def _bits(values: np.ndarray) -> list[list[int]]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def _inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    image = nib.load(str(ROOT / "optimized_centroid_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 4, 3, 12)
    return (
        frames[..., 0].astype(np.int16),
        frames[..., 1:6],
        frames[..., 6:8],
        frames[..., 8:10],
    )


def test_optimized_tissue_centroid_modes_match_java_powertable() -> None:
    assert REFERENCE["source_sha256"] == "77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0"
    assert REFERENCE["powertable_sha256"] == POWER_TABLE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    labels, memberships, images, fields = _inputs()
    assert set(np.unique(labels)) == {0, 1, 2, 3}  # fifth class is absent
    # c0 at (1, 0, 0) has an exact zero residual even after its unit field;
    # c0 at (0, 0, 2) is a no-field residual below the source ZERO cutoff.
    assert labels[1, 0, 0] == 1 and fields[1, 0, 0, 0] == np.float32(1.0)
    assert images[1, 0, 0, 0] == INITIAL[0, 1]
    assert labels[0, 0, 2] == 2
    tiny_residual = np.float32(images[0, 0, 2, 0] - INITIAL[0, 2])
    assert np.float32(tiny_residual * tiny_residual) < np.float32(1e-10)
    for mode in ("none", "linked", "damped", "prior", "fallback"):
        for field_key, field_values in (("no_field", None), ("with_field", fields)):
            kwargs: dict[str, object] = {
                "fields": field_values,
                "mode": mode,
                "smoothness": 0.35,
                "power_exponents": EXPONENTS,
                "modalities": MODALITIES,
            }
            if mode == "linked":
                kwargs["centroid_links"] = LINKS
            if mode == "prior":
                kwargs["prior_centroids"] = PRIOR
            actual = compute_lesion_toads_basic_centroids(
                images, memberships, labels, INITIAL, TYPES, **kwargs
            )
            assert _bits(actual) == REFERENCE[f"{mode}_{field_key}"]


def test_optimized_centroids_reject_bad_exponents_and_modalities() -> None:
    labels, memberships, images, _fields = _inputs()
    with pytest.raises(ValueError, match="power_exponents"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES,
            power_exponents=np.zeros((2, 4), dtype=np.float32), modalities=MODALITIES,
        )
    with pytest.raises(ValueError, match="modalities"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES,
            power_exponents=EXPONENTS, modalities=("T1_MPRAGE",),
        )
    with pytest.raises(ValueError, match="finite"):
        invalid = EXPONENTS.copy()
        invalid[0, 1] = np.nan
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES,
            power_exponents=invalid, modalities=MODALITIES,
        )
    with pytest.raises(ValueError, match="nonnegative"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES,
            power_exponents=EXPONENTS, modalities=MODALITIES, zero_threshold=-1,
        )
