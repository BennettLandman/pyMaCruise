from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_data_energy_maps import compute_lesion_toads_data_energy_maps


FIXTURE = Path(__file__).parent / "data/synthetic/lesion_data_energy_maps"
HARNESS = Path(__file__).parent / "java_reference/LesionDataEnergyMapsReferenceHarness.java"
IMAGE_SHA256 = "cc7a4e4ccfea1094f7c09475d62f204d0a972c865a48dc763d38d42f9435d0fc"
REFERENCE_SHA256 = "658b2f1df7dc4c5f0125596d7025906024d91546fd2010ed438ebee744c9efb4"
HARNESS_SHA256 = "82a4db72e6db53280f19d7ec5f80df4c02c01c4bb4faf41c3e806f0a8f449853"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA256 = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"


def _from_bits(values: list[int]) -> np.ndarray:
    return np.asarray([value & 0xFFFFFFFF for value in values], dtype=np.uint32).view(np.float32)


def _bits(values: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=np.float32).ravel().view(np.uint32)


def _fixture():
    path = FIXTURE / "images_ras.nii.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == IMAGE_SHA256
    image = nib.load(path)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert np.array_equal(nib.affines.voxel_sizes(image.affine), (1.5, 2.0, 2.5))
    values = np.asarray(image.dataobj, dtype=np.float32)
    assert values.shape == (5, 6, 7, 3)
    return values


def _inputs(case: dict, values: np.ndarray):
    channels = values.shape[-1]
    classes = len(case["class_types"])
    centroids = _from_bits(case["class_centroid_bits"]).reshape(classes, channels)
    lesion = _from_bits(case["lesion_centroid_bits"])
    black_hole = _from_bits(case["black_hole_centroid_bits"])
    weights = _from_bits(case["channel_weight_bits"])
    return centroids, lesion, black_hole, weights


def test_volume_energy_maps_match_selected_java_and_ras_fixture() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    reference_path = FIXTURE / "java_reference.json"
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == REFERENCE_SHA256
    reference = json.loads(reference_path.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["wrapper_sha256"] == WRAPPER_SHA256
    assert reference["shape"] == [5, 6, 7]
    assert reference["channels"] == 3

    values = _fixture()
    # The RAS NIfTI samples match the deterministic values in the standalone
    # Java harness, including each channel at asymmetric spatial coordinates.
    x, y, z, c = np.indices(values.shape, dtype=np.int32)
    source_values = (13 * x - 7 * y + 5 * z + 11 * c - 23).astype(np.float32) / np.float32(8.0)
    np.testing.assert_array_equal(values, source_values)
    original = values.copy()

    for case in reference["cases"]:
        centroids, lesion, black_hole, weights = _inputs(case, values)
        result = compute_lesion_toads_data_energy_maps(
            values,
            centroids,
            lesion,
            black_hole,
            weights,
            case["class_types"],
        )
        shape = tuple(reference["shape"])
        nvox = int(np.prod(shape))
        nclasses = len(case["class_types"])
        assert result.class_energy.shape == (*shape, nclasses)
        assert result.lesion_energy.shape == shape
        assert result.black_hole_energy.shape == shape
        np.testing.assert_array_equal(
            _bits(result.class_energy),
            np.asarray([v & 0xFFFFFFFF for v in case["class_energy_bits"]], dtype=np.uint32),
            err_msg=case["name"],
        )
        np.testing.assert_array_equal(
            _bits(result.lesion_energy),
            np.asarray([v & 0xFFFFFFFF for v in case["lesion_energy_bits"]], dtype=np.uint32),
            err_msg=case["name"],
        )
        np.testing.assert_array_equal(
            _bits(result.black_hole_energy),
            np.asarray([v & 0xFFFFFFFF for v in case["black_hole_energy_bits"]], dtype=np.uint32),
            err_msg=case["name"],
        )

        # Pin the overwrite rule explicitly: later auxiliary classes win the
        # shared lesion slot; a later LESION leaves the prior BH slot intact.
        lesion_by_class = np.asarray(
            [v & 0xFFFFFFFF for v in case["lesion_energy_by_class_bits"]], dtype=np.uint32
        ).reshape(nvox, nclasses)
        black_hole_by_class = np.asarray(
            [v & 0xFFFFFFFF for v in case["black_hole_energy_by_class_bits"]], dtype=np.uint32
        ).reshape(nvox, nclasses)
        last_lesion_related = max(
            i for i, kind in enumerate(case["class_types"]) if kind in {"LESION", "BLACKHOLE"}
        )
        last_black_hole = max(i for i, kind in enumerate(case["class_types"]) if kind == "BLACKHOLE")
        result_lesion = _bits(result.lesion_energy)
        result_bh = _bits(result.black_hole_energy)
        # Compare interior positions using x/y/z C-order indexing.
        for ix in range(1, shape[0] - 1):
            for iy in range(1, shape[1] - 1):
                for iz in range(1, shape[2] - 1):
                    flat_voxel = (ix * shape[1] + iy) * shape[2] + iz
                    assert result_lesion[flat_voxel] == lesion_by_class[flat_voxel, last_lesion_related]
                    assert result_bh[flat_voxel] == black_hole_by_class[flat_voxel, last_black_hole]
        np.testing.assert_array_equal(values, original)

        boundary = np.ones(shape, dtype=bool)
        boundary[1:-1, 1:-1, 1:-1] = False
        np.testing.assert_array_equal(result.class_energy[boundary], 0.0)
        np.testing.assert_array_equal(result.lesion_energy[boundary], 0.0)
        np.testing.assert_array_equal(result.black_hole_energy[boundary], 0.0)


def test_zero_residual_energy_maps_are_zero_on_empty_signal() -> None:
    images = np.zeros((3, 3, 3, 2), dtype=np.float32)
    result = compute_lesion_toads_data_energy_maps(
        images,
        np.zeros((3, 2), dtype=np.float32),
        np.zeros(2, dtype=np.float32),
        np.zeros(2, dtype=np.float32),
        np.asarray([0.5, 1.5], dtype=np.float32),
        ["OBJECT", "LESION", "BLACKHOLE"],
    )
    assert not np.any(result.class_energy)
    assert not np.any(result.lesion_energy)
    assert not np.any(result.black_hole_energy)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"class_types": ["OBJECT"]},
        {"class_types": ["OBJECT", "RAYLEIGH"]},
    ],
)
def test_volume_energy_maps_reject_unsupported_class_layout(kwargs: dict) -> None:
    images = np.zeros((3, 3, 3, 1), dtype=np.float32)
    with pytest.raises(ValueError):
        compute_lesion_toads_data_energy_maps(
            images,
            np.zeros((2, 1), dtype=np.float32),
            np.zeros(1, dtype=np.float32),
            np.zeros(1, dtype=np.float32),
            np.ones(1, dtype=np.float32),
            kwargs["class_types"],
        )
