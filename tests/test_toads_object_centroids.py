from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_object_centroids import toads_object_centroid_update


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_object_centroids"
JAVA = Path(__file__).parent / "java_reference" / "ToadsObjectCentroidsReferenceHarness.java"
REFERENCE = json.loads((Path(__file__).parent / "java_reference" / "ToadsObjectCentroidsReference.json").read_text())
HARNESS_SHA256 = "4532710f49ff5455336de69b88ee0bf9972b3adf66c276fcdea8164ffcfb497a"
WRAPPER_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
INITIAL = np.asarray([[0.31, 0.72, 1.19], [0.43, 0.88, 1.37]], dtype=np.float32)


def _inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    frames_image = nib.load(str(ROOT / "object_centroid_inputs.nii.gz"))
    segmentation_image = nib.load(str(ROOT / "segmentation.nii.gz"))
    assert nib.aff2axcodes(frames_image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(segmentation_image.affine) == ("R", "A", "S")
    values = np.asarray(frames_image.dataobj, dtype=np.float32)
    return values[..., :2], values[..., 2:4], values[..., 4:7], np.asarray(segmentation_image.dataobj)


def _bits(values: np.ndarray) -> list[list[int]]:
    return values.astype(np.float32, copy=False).view(np.uint32).astype(np.uint64).tolist()


def test_ordinary_object_centroids_match_selected_java_excerpt_bits() -> None:
    assert hashlib.sha256(JAVA.read_bytes()).hexdigest() == HARNESS_SHA256
    assert REFERENCE["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert REFERENCE["wrapper_sha256"] == WRAPPER_SHA256
    images, field, memberships, segmentation = _inputs()
    assert images.shape == field.shape == (5, 4, 4, 2)
    assert memberships.shape == (5, 4, 4, 3)
    assert np.all(memberships[..., 2] == 0)

    cases = {
        "raw": toads_object_centroid_update(images, memberships, INITIAL),
        "field": toads_object_centroid_update(
            images, memberships, INITIAL, intensity_field=field
        ),
        "ridge": toads_object_centroid_update(
            images, memberships, INITIAL, segmentation=segmentation, psf_ridge=0.23
        ),
        "field_ridge": toads_object_centroid_update(
            images,
            memberships,
            INITIAL,
            intensity_field=field,
            segmentation=segmentation,
            psf_ridge=0.23,
        ),
    }
    for case, actual in cases.items():
        assert _bits(actual) == REFERENCE[case]
        assert _bits(actual[:, 2:3]) == _bits(INITIAL[:, 2:3])


def test_centroid_update_requires_halo_and_segmentation_for_ridge() -> None:
    images, _, memberships, segmentation = _inputs()
    with pytest.raises(ValueError, match="one-voxel halo"):
        toads_object_centroid_update(images[:2], memberships[:2], INITIAL)
    with pytest.raises(ValueError, match="segmentation"):
        toads_object_centroid_update(images, memberships, INITIAL, psf_ridge=0.23)
    with pytest.raises(ValueError, match="smoothness"):
        toads_object_centroid_update(images, memberships, INITIAL, smoothness=1.1)
    with pytest.raises(ValueError, match="intensity_field"):
        toads_object_centroid_update(
            images, memberships, INITIAL, intensity_field=segmentation[..., None]
        )
