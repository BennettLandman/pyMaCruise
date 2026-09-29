from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_shapes import sample_lesion_shape_prior


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_shape_sampling"
FIXTURE = FIXTURE_DIR / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "LesionShapeSamplingReferenceHarness.java"
HARNESS_SHA256 = "b5e32e1e77eb66c5edca71b275160963d9d1590d2b77d3e5ccba6a550c78e075"


def _shape_prior(dimensions: tuple[int, int, int], class_index: int) -> np.ndarray:
    offset = np.float32(0.12 if class_index == 0 else 0.23)
    shape = np.empty(dimensions, dtype=np.float32)
    for x in range(dimensions[0]):
        for y in range(dimensions[1]):
            for z in range(dimensions[2]):
                value = np.float32(offset + np.float32(np.float32(0.017) * x))
                value = np.float32(value + np.float32(np.float32(0.023) * y))
                value = np.float32(value + np.float32(np.float32(0.031) * z))
                value = np.float32(
                    value + np.float32(np.float32(np.float32(0.004) * x) * y)
                )
                shape[x, y, z] = value
    return shape


def test_lesion_shape_sampler_matches_selected_java_formula_harness() -> None:
    reference = json.loads(FIXTURE.read_text())
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for case in reference["cases"]:
        sample = sample_lesion_shape_prior(
            _shape_prior(tuple(reference["template_dims"]), case["class"]),
            tuple(case["point"]),
            image_dimensions_xyz=tuple(reference["image_dims"]),
            image_resolutions_mm=tuple(reference["image_resolutions"]),
            template_resolutions_mm=tuple(reference["template_resolutions"]),
            transform=tuple(case["transform"]),
            class_index=case["class"],
            bbox_low_xyz=tuple(case["bbox_low"]),
            bbox_high_xyz=tuple(case["bbox_high"]),
        )
        assert sample.value == np.float32(case["value"]), case["name"]
        if case["mapped"] is None:
            assert sample.template_xyz is None
        else:
            np.testing.assert_array_equal(
                np.asarray(sample.template_xyz, dtype=np.float32),
                np.asarray(case["mapped"], dtype=np.float32),
                err_msg=case["name"],
            )


def test_lesion_shape_none_mode_directly_samples_ras_loaded_nifti_array(tmp_path: Path) -> None:
    values = np.arange(4 * 5 * 6, dtype=np.float32).reshape((4, 5, 6)) / 100.0
    affine = np.diag([0.8, 1.2, 1.7, 1.0])
    path = tmp_path / "shape.nii.gz"
    nib.save(nib.Nifti1Image(values, affine), path)
    loaded = nib.load(str(path))
    sample = sample_lesion_shape_prior(
        np.asarray(loaded.dataobj, dtype=np.float32),
        (2, 3, 4),
        image_dimensions_xyz=(4, 5, 6),
        image_resolutions_mm=(0.8, 1.2, 1.7),
        template_resolutions_mm=(0.8, 1.2, 1.7),
        transform=(0, 0, 0, 0, 0, 0),
        class_index=2,
        bbox_low_xyz=(0, 0, 0),
        bbox_high_xyz=(4, 5, 6),
        transform_mode="none",
    )
    assert nib.aff2axcodes(loaded.affine) == ("R", "A", "S")
    assert sample.value == values[2, 3, 4]
    assert sample.template_xyz is None


def test_lesion_shape_sampler_rejects_unsupported_modes_and_bad_maps() -> None:
    arguments = {
        "image_dimensions_xyz": (4, 4, 4),
        "image_resolutions_mm": (1.0, 1.0, 1.0),
        "template_resolutions_mm": (1.0, 1.0, 1.0),
        "transform": (0, 0, 0, 0, 0, 0),
        "class_index": 1,
        "bbox_low_xyz": (0, 0, 0),
        "bbox_high_xyz": (4, 4, 4),
    }
    with pytest.raises(ValueError, match="only 'none' or 'rigid'"):
        sample_lesion_shape_prior(
            np.ones((4, 4, 4), dtype=np.float32),
            (1, 1, 1),
            transform_mode="deformable",
            **arguments,
        )
    with pytest.raises(ValueError, match="match the input grid"):
        sample_lesion_shape_prior(
            np.ones((3, 3, 3), dtype=np.float32),
            (1, 1, 1),
            transform_mode="none",
            **arguments,
        )
    with pytest.raises(ValueError, match="finite"):
        sample_lesion_shape_prior(
            np.full((4, 4, 4), np.nan, dtype=np.float32),
            (1, 1, 1),
            **arguments,
        )
