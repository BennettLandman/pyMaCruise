from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_multichannel_exact_field_refit import toads_multichannel_exact_field_refit
from test_toads_post_thinning_exact_iteration import FIELD_HARNESS, _field_input_bytes

ROOT = Path(__file__).parent
CASE_DIR = ROOT / "data/synthetic/toads_p18czo_multichannel_exact_field"
CASE = json.loads((CASE_DIR / "case.json").read_text())
SOURCE_DIR = ROOT / "data/synthetic" / CASE["source_fixture_directory"]
JAVA_FIELD_SOURCE_SHA256 = "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8"


def _load(path: Path) -> nib.Nifti1Image:
    return nib.load(str(path))


def test_two_channel_degree_one_field_refit_matches_selected_java_formulas(tmp_path: Path) -> None:
    for name, expected_sha in CASE["fixture_sha256"].items():
        path = CASE_DIR / name if name == CASE["segmentation"] else SOURCE_DIR / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_sha
    image_volumes = (
        _load(SOURCE_DIR / CASE["channel_order"][0].split(" / ")[0]),
        _load(SOURCE_DIR / CASE["channel_order"][1].split(" / ")[0]),
    )
    field_volumes = (
        _load(SOURCE_DIR / CASE["channel_order"][0].split(" / ")[1]),
        _load(SOURCE_DIR / CASE["channel_order"][1].split(" / ")[1]),
    )
    membership_volume = _load(SOURCE_DIR / "exact_memberships.nii.gz")
    segmentation_volume = _load(CASE_DIR / CASE["segmentation"])
    all_volumes = (*image_volumes, *field_volumes, membership_volume, segmentation_volume)
    reference_affine = image_volumes[0].affine
    for volume in all_volumes:
        assert nib.aff2axcodes(volume.affine) == tuple(CASE["axis_codes"])
        np.testing.assert_allclose(volume.affine, reference_affine, rtol=0, atol=1e-6)
    np.testing.assert_allclose(nib.affines.voxel_sizes(reference_affine), CASE["voxel_spacing_mm"], atol=1e-6)

    images = np.stack(tuple(np.asarray(volume.dataobj, dtype=np.float32) for volume in image_volumes), axis=-1)
    fields = np.stack(tuple(np.asarray(volume.dataobj, dtype=np.float32) for volume in field_volumes), axis=-1)
    memberships = np.asarray(membership_volume.dataobj, dtype=np.float32)
    segmentation = np.asarray(segmentation_volume.dataobj, dtype=np.int16)
    centers = np.stack((
        np.asarray([1.0 / 12.0 + k / 12.0 for k in range(11)], dtype=np.float32),
        np.asarray([(11 - k) / 12.0 for k in range(11)], dtype=np.float32),
    ))
    original_fields = fields.copy()
    result = toads_multichannel_exact_field_refit(
        images, memberships, segmentation, centers, fields
    )
    assert result.channel_order == (0, 1)
    assert len(result.fits) == CASE["channels"] == 2
    assert result.fields.shape == images.shape
    np.testing.assert_array_equal(fields, original_fields)

    for channel in range(CASE["channels"]):
        input_path = tmp_path / f"channel_{channel}_field_fit.bin"
        input_path.write_bytes(_field_input_bytes(
            segmentation,
            images[..., channel],
            memberships,
            original_fields[..., channel],
            centers[channel],
        ))
        reference = json.loads(subprocess.run(
            ["sh", str(FIELD_HARNESS), str(input_path), *(str(v) for v in segmentation.shape), "11"],
            check=True, capture_output=True, text=True,
        ).stdout)
        assert reference["reference_kind"] == "formula_harness_not_source_class_execution"
        assert reference["inhomogeneity_upstream_source_sha256"] == JAVA_FIELD_SOURCE_SHA256
        java_field = np.asarray(reference["field_bits"], dtype=np.int32).view(np.float32).reshape(segmentation.shape)
        np.testing.assert_allclose(result.fields[..., channel], java_field, rtol=0, atol=2e-7)
        np.testing.assert_allclose(result.fits[channel].coefficients, reference["coefficients"], rtol=2e-6, atol=2e-8)
        np.testing.assert_allclose(result.fits[channel].normal_matrix, reference["normal"], rtol=2e-6, atol=2e-8)
        np.testing.assert_allclose(result.fits[channel].right_hand_side, reference["rhs"], rtol=2e-6, atol=2e-8)
        assert not np.array_equal(result.fields[..., channel], original_fields[..., channel])

    assert not np.array_equal(result.fields[..., 0], result.fields[..., 1])
