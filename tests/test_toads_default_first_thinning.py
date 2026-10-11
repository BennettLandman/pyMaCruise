from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_first_thinning import toads_first_exact_thinning_pass
from pymacruise.toads_thinning import toads_initialize_thinning_queue


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_first_thinning"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _float_digest(values: np.ndarray) -> str:
    return _sha256(np.asarray(values, dtype=np.float32).view(np.uint32).astype(">u4").tobytes())


def _label_digest(values: np.ndarray) -> str:
    return _sha256(np.asarray(values, dtype=np.int8).tobytes())


def _queue_digest(items) -> str:
    trace = "".join(
        f"{int(np.asarray([item.value], dtype=np.float32).view(np.int32)[0])},"
        f"{item.x},{item.y},{item.z},{item.label}\n"
        for item in items
    )
    return _sha256(trace.encode("ascii"))


def _available_digest(values: np.ndarray) -> str:
    text = "".join("1" if value else "0" for value in np.asarray(values, dtype=bool).ravel())
    return _sha256(text.encode("ascii"))


def test_selected_exact_initial_memberships_feed_first_default_thinning_pass() -> None:
    spec = json.loads((ROOT / "spec.json").read_text())
    image = nib.load(str(ROOT / "subject_images.nii.gz"))
    segmentation_image = nib.load(str(ROOT / "initial_segmentation.nii.gz"))
    assert image.shape == tuple(spec["shape_xyz_channels"])
    assert segmentation_image.shape == tuple(spec["shape_xyz_channels"][:3])
    assert nib.aff2axcodes(image.affine) == nib.aff2axcodes(segmentation_image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(image.affine, segmentation_image.affine)
    assert image.header.get_xyzt_units()[0] == "mm"

    images = np.moveaxis(np.asarray(image.dataobj, dtype=np.float32), -1, 0)
    segmentation = np.asarray(segmentation_image.dataobj, dtype=np.int16)
    centroids = np.asarray(spec["centroids"], dtype=np.float32)
    classes = spec["classes"]
    relations = np.zeros((classes, classes), dtype=bool)
    shapes = np.zeros((*segmentation.shape, classes), dtype=np.float32)

    # This is the explicit selected default ctor profile: full exact class
    # memberships, not the alternate approx=4 constructor.
    result = toads_first_exact_thinning_pass(
        images,
        centroids,
        segmentation,
        relations,
        shapes,
        crop_border=spec["crop_border"],
        volume_ratio=spec["volume_ratio"],
        limit=spec["thinning_limit"],
        # This selected Java comparison harness explicitly initializes the
        # skeleton to the empty-label sentinel before initThinningLabels.
        initial_skeleton=np.full(segmentation.shape, -1, dtype=np.int16),
    )
    initial = toads_initialize_thinning_queue(
        segmentation,
        result.memberships,
        zero_threshold=1e-10,
    )
    thinning = result.thinning

    assert REFERENCE["source"] == (
        "VandyCRUISE MedicAlgorithmToads08 exact initial memberships -> "
        "MultipleToadSegmentation.initThinningLabels + propagateThinning"
    )
    assert REFERENCE["preprocess_source_sha256"] == "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7"
    assert REFERENCE["segmentation_source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert REFERENCE["algorithm_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert REFERENCE["wrapper_source_sha256"] == "d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115"
    assert REFERENCE["lut_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert REFERENCE["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert REFERENCE["shape"] == [36, 36, 36, 11]
    assert REFERENCE["crop_border"] == 10

    assert _float_digest(result.memberships) == REFERENCE["memberships_sha256"]
    np.testing.assert_array_equal(_bits(initial.class_volumes), REFERENCE["initial_volume_bits"])
    np.testing.assert_array_equal(_bits(thinning.minimum_class_volumes), REFERENCE["minimum_volume_bits"])
    assert initial.boundary_counts.tolist() == REFERENCE["boundary_counts"]
    assert len(thinning.processed_items) == REFERENCE["processed_count"]
    assert _queue_digest(thinning.processed_items) == REFERENCE["processed_sha256"]
    assert _label_digest(thinning.segmentation) == REFERENCE["segmentation_sha256"]
    assert _label_digest(thinning.skeleton) == REFERENCE["skeleton_sha256"]
    # This source-selected wrapper pre-seeds the full skeleton with -1, so
    # the crop border remains empty after initThinningLabels and propagation.
    assert np.all(thinning.skeleton[0, :, :] == -1)
    assert np.all(thinning.skeleton[-1, :, :] == -1)
    assert np.all(thinning.skeleton[:, 0, :] == -1)
    assert np.all(thinning.skeleton[:, -1, :] == -1)
    assert np.all(thinning.skeleton[:, :, 0] == -1)
    assert np.all(thinning.skeleton[:, :, -1] == -1)
    assert _available_digest(thinning.available) == REFERENCE["available_sha256"]
    np.testing.assert_array_equal(_bits(thinning.class_volumes), REFERENCE["final_volume_bits"])
    assert len(thinning.queue) == REFERENCE["queue_after_count"]
    assert REFERENCE["queue_after_sha256"] == _sha256(b"")

    # All selected-class volumes exceed the source 100-voxel floor before
    # thinning; default wrapper vLim=.5 supplies the larger fractional floor.
    assert np.all(initial.class_volumes > np.float32(100.0))
    assert np.all(thinning.minimum_class_volumes > np.float32(100.0))
    assert np.any(thinning.class_volumes < thinning.minimum_class_volumes)
