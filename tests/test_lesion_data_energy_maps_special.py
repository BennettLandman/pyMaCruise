from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_data_energy_maps import compute_lesion_toads_data_energy_maps


FIXTURE = Path(__file__).parent / "data/synthetic/lesion_data_energy_maps_special"
HARNESS = Path(__file__).parent / "java_reference/LesionDataEnergyMapsSpecialReferenceHarness.java"
IMAGE_SHA256 = "6fece9958c7040cc65c8f00e102ad97e0bd997a0e6a3705f07ee39d2b3b2e2e2"
REFERENCE_SHA256 = "c451338a135ca8b3d87a0f6106e24f784aa95df32439f2e241013961b351aad4"
HARNESS_SHA256 = "82bda341b858702507d58aa236327d37c08b8a95302ffa6445590ca9f218101c"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"


def _from_bits(values: list[int] | int) -> np.ndarray | np.float32:
    if isinstance(values, int):
        return np.asarray(values & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]
    return np.asarray([value & 0xFFFFFFFF for value in values], dtype=np.uint32).view(np.float32)


def _bits(values: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=np.float32).ravel().view(np.uint32)


def _load_inputs():
    path = FIXTURE / "inputs_ras.nii.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == IMAGE_SHA256
    nii = nib.load(path)
    assert nib.aff2axcodes(nii.affine) == ("R", "A", "S")
    assert np.array_equal(nib.affines.voxel_sizes(nii.affine), (1.5, 2.0, 2.5))
    data = np.asarray(nii.dataobj, dtype=np.float32)
    assert data.shape == (7, 8, 9, 6)
    return data[..., :3], data[..., 3:]


def _check_fixture_formulas(images: np.ndarray, fields: np.ndarray) -> None:
    x, y, z, c = np.indices(images.shape, dtype=np.int32)
    fold = (3 * x + 2 * y + z) % 5 - 2
    expected_image = (13 * x - 9 * y + 7 * z + 5 * fold + 11 * c - 20).astype(np.float32) / np.float32(32.0)
    expected_field = np.float32(0.625) + ((x + 2 * y + 3 * z + c) % 7).astype(np.float32) / np.float32(16.0)
    np.testing.assert_array_equal(images, expected_image)
    np.testing.assert_array_equal(fields, expected_field)


def test_fielded_special_volume_maps_match_selected_java_harness() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    reference_path = FIXTURE / "java_reference.json"
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == REFERENCE_SHA256
    reference = json.loads(reference_path.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["power_table_sha256"] == POWER_TABLE_SHA256
    assert reference["shape"] == [7, 8, 9]
    assert reference["channels"] == 3
    assert len(reference["cases"]) == 9

    images, fields = _load_inputs()
    _check_fixture_formulas(images, fields)
    originals = images.copy(), fields.copy()
    outputs = {}
    for case in reference["cases"]:
        channels = images.shape[-1]
        nclasses = len(case["class_types"])
        result = compute_lesion_toads_data_energy_maps(
            images,
            _from_bits(case["class_centroid_bits"]).reshape(nclasses, channels),
            _from_bits(case["lesion_centroid_bits"]),
            _from_bits(case["black_hole_centroid_bits"]),
            _from_bits(case["channel_weight_bits"]),
            case["class_types"],
            use_field=case["use_field"],
            intensity_fields=fields,
            structure_names=case["structure_names"],
            t1_channel=case["t1_channel"],
            has_mprage=case["has_mprage"],
            csf_power_exponent=float(_from_bits(case["csf_power_exponent_bits"])),
            flair_channel=case["flair_channel"],
            use_lesion_weight=case["use_lesion_weight"],
            lesion_power_exponent=float(_from_bits(case["lesion_power_exponent_bits"])),
        )
        shape = tuple(reference["shape"])
        nvox = int(np.prod(shape))
        assert result.class_energy.shape == (*shape, nclasses)
        np.testing.assert_array_equal(
            _bits(result.class_energy),
            np.asarray([value & 0xFFFFFFFF for value in case["class_energy_bits"]], dtype=np.uint32),
            err_msg=case["name"],
        )
        np.testing.assert_array_equal(
            _bits(result.lesion_energy),
            np.asarray([value & 0xFFFFFFFF for value in case["lesion_energy_bits"]], dtype=np.uint32),
            err_msg=case["name"],
        )
        np.testing.assert_array_equal(
            _bits(result.black_hole_energy),
            np.asarray([value & 0xFFFFFFFF for value in case["black_hole_energy_bits"]], dtype=np.uint32),
            err_msg=case["name"],
        )

        # The BH class precedes LESION: lesion energy is replaced by the
        # later LESION FLAIR branch while the BH slot retains its earlier value.
        lesion_by_class = np.asarray(
            [value & 0xFFFFFFFF for value in case["lesion_energy_by_class_bits"]], dtype=np.uint32
        ).reshape(nvox, nclasses)
        bh_by_class = np.asarray(
            [value & 0xFFFFFFFF for value in case["black_hole_energy_by_class_bits"]], dtype=np.uint32
        ).reshape(nvox, nclasses)
        lesion_result = _bits(result.lesion_energy)
        bh_result = _bits(result.black_hole_energy)
        for x in range(1, shape[0] - 1):
            for y in range(1, shape[1] - 1):
                for z in range(1, shape[2] - 1):
                    voxel = (x * shape[1] + y) * shape[2] + z
                    assert lesion_result[voxel] == lesion_by_class[voxel, 3]
                    assert bh_result[voxel] == bh_by_class[voxel, 2]

        boundary = np.ones(shape, dtype=bool)
        boundary[1:-1, 1:-1, 1:-1] = False
        np.testing.assert_array_equal(result.class_energy[boundary], 0.0)
        np.testing.assert_array_equal(result.lesion_energy[boundary], 0.0)
        np.testing.assert_array_equal(result.black_hole_energy[boundary], 0.0)
        outputs[case["name"]] = result

    np.testing.assert_array_equal(images, originals[0])
    np.testing.assert_array_equal(fields, originals[1])
    assert not np.array_equal(
        outputs["unfielded_sulcal_direct"].class_energy,
        outputs["fielded_sulcal_direct"].class_energy,
    )
    assert not np.array_equal(
        outputs["unfielded_sulcal_direct"].class_energy,
        outputs["unfielded_sulcal_mprage"].class_energy,
    )
    assert not np.array_equal(
        outputs["lesion_flair_unweighted"].lesion_energy,
        outputs["lesion_flair_power_unfielded"].lesion_energy,
    )


def test_volume_field_and_channel_options_are_validated() -> None:
    images = np.zeros((3, 3, 3, 2), dtype=np.float32)
    common = (
        images,
        np.zeros((1, 2), dtype=np.float32),
        np.zeros(2, dtype=np.float32),
        np.zeros(2, dtype=np.float32),
        np.ones(2, dtype=np.float32),
        ["OBJECT"],
    )
    with pytest.raises(ValueError, match="intensity_fields are required"):
        compute_lesion_toads_data_energy_maps(*common, use_field=True)
    with pytest.raises(ValueError, match="t1_channel"):
        compute_lesion_toads_data_energy_maps(*common, t1_channel=2)
    with pytest.raises(ValueError, match="structure_names"):
        compute_lesion_toads_data_energy_maps(*common, structure_names=["OBJECT", "extra"])
