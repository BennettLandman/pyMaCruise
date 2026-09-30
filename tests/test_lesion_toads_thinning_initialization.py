from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import initialize_lesion_toads_thinning_frontier


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_thinning_initialization"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
OBJECT_TYPES = ("obj", "lesion", "blackhole")


def _load(name: str, dtype=None) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj, dtype=dtype)


def _u32(values) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _reference_u32(values) -> np.ndarray:
    return np.asarray(values, dtype=np.int64) & 0xFFFFFFFF


def _input_arrays():
    return (
        _load("segmentation", np.int16),
        _load("memberships", np.float32),
        _load("lesion_membership", np.float32),
        _load("black_hole_membership", np.float32),
        _load("initial_skeleton", np.int16),
        _load("initial_available").astype(bool),
    )


def test_thinning_frontier_initialization_matches_selected_java_and_heap() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert REFERENCE["heap_source_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert REFERENCE["shape"] == [4, 4, 4]
    arrays = _input_arrays()
    seg, memberships, lesion, black_hole, skeleton, available = arrays
    original_arrays = tuple(array.copy() for array in arrays)

    x, y, z = np.indices(seg.shape)
    expected_segmentation = ((x + y + z) % 3).astype(np.int16)
    for point in ((0, 1, 1), (2, 1, 1), (1, 0, 1), (1, 2, 1), (1, 1, 0), (1, 1, 2)):
        expected_segmentation[point] = 0
    np.testing.assert_array_equal(seg, expected_segmentation)
    expected_memberships = np.stack(
        [
            (2 + ((x * 7 + y * 5 + z * 3 + k * 11) % 15)).astype(np.float32)
            / np.float32(32.0)
            for k in range(3)
        ],
        axis=-1,
    )
    np.testing.assert_array_equal(memberships, expected_memberships)

    for reference in REFERENCE["cases"]:
        result = initialize_lesion_toads_thinning_frontier(
            seg,
            memberships,
            OBJECT_TYPES,
            lesion,
            black_hole,
            skeleton,
            available,
            lesion_class=1,
            use_lesion_max=reference["use_lesion_max"],
        )
        np.testing.assert_array_equal(
            _u32(result.class_volume), _reference_u32(reference["volume_bits"])
        )
        np.testing.assert_array_equal(result.boundary_counts, reference["boundary_counts"])
        np.testing.assert_array_equal(result.skeleton.ravel(), reference["skeleton"])
        np.testing.assert_array_equal(
            result.available.ravel().tolist(), reference["available"]
        )
        np.testing.assert_array_equal(
            _u32([item.value for item in result.inserted_items]),
            _reference_u32(reference["inserted_priority_bits"]),
        )
        np.testing.assert_array_equal(
            np.asarray(
                [[item.x, item.y, item.z, item.label] for item in result.inserted_items],
                dtype=np.int64,
            ).ravel(),
            reference["inserted_payloads"],
        )
        popped = [result.queue.pop_first() for _ in range(result.queue.size)]
        np.testing.assert_array_equal(
            _u32([item.value for item in popped]),
            _reference_u32(reference["removed_priority_bits"]),
        )
        np.testing.assert_array_equal(
            np.asarray(
                [[item.x, item.y, item.z, item.label] for item in popped], dtype=np.int64
            ).ravel(),
            reference["removed_payloads"],
        )

    for actual, original in zip(arrays, original_arrays, strict=True):
        np.testing.assert_array_equal(actual, original)


def test_thinning_frontier_initializer_rejects_class_grid_mismatch() -> None:
    seg, memberships, lesion, black_hole, skeleton, available = _input_arrays()
    with np.testing.assert_raises_regex(ValueError, "boolean grid"):
        initialize_lesion_toads_thinning_frontier(
            seg,
            memberships,
            OBJECT_TYPES,
            lesion,
            black_hole,
            skeleton,
            available[..., :2],
            lesion_class=1,
        )
