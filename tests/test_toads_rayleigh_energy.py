from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_rayleigh_energy import toads_rayleigh_intensity_energy


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_rayleigh_energy"
SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _from_bits(bits: list[int]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32)


def _bits(values: np.ndarray) -> list[int]:
    canonical = np.asarray(values, dtype=np.float32).copy()
    canonical.view(np.uint32)[np.isnan(canonical)] = np.uint32(0x7FC00000)
    return canonical.view(np.int32).ravel().tolist()


def _fixture() -> tuple[dict, dict]:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    image = nib.load(str(ROOT / "images.nii.gz"))
    field = nib.load(str(ROOT / "field.nii.gz"))
    bias = nib.load(str(ROOT / "psf_bias.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), [0.85, 1.15, 1.4], rtol=0, atol=1e-6)
    np.testing.assert_allclose(field.affine, image.affine, rtol=0, atol=1e-7)
    np.testing.assert_allclose(bias.affine, image.affine, rtol=0, atol=1e-7)
    images = np.asarray(image.dataobj, dtype=np.float32)
    fields = np.asarray(field.dataobj, dtype=np.float32)
    biases = np.asarray(bias.dataobj, dtype=np.float32)
    assert images.shape == fields.shape == biases.shape == (11, 9, 7, 3)
    params = {
        "centroids": _from_bits(reference["centroid_bits"]).reshape(3, 4),
        "modality_weights": _from_bits(reference["modality_bits"]),
        "class_variance_weights": _from_bits(reference["class_variance_bits"]).reshape(3, 4),
        "image_variances": _from_bits(reference["image_variance_bits"]),
    }
    return {"image_values": images, "field": fields, "bias": biases}, params


def test_rayleigh_class_energy_matches_selected_java_all_modes() -> None:
    volumes, params = _fixture()
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["shape"] == [11, 9, 7, 4]
    inputs = {name: value.copy() for name, value in volumes.items()}
    for case in reference["cases"]:
        actual = toads_rayleigh_intensity_energy(
            volumes["image_values"],
            **params,
            intensity_field=volumes["field"] if case["field"] else None,
            psf_bias=volumes["bias"] if case["bias"] else None,
        )
        assert actual.shape == (11, 9, 7, 4)
        assert actual.dtype == np.float32
        assert _bits(actual) == case["energy_bits"]
        # The adapter materializes only the source's visited interior; all
        # six faces remain at their zero-initialized values.
        assert np.all(actual[0] == 0) and np.all(actual[-1] == 0)
        assert np.all(actual[:, 0] == 0) and np.all(actual[:, -1] == 0)
        assert np.all(actual[:, :, 0] == 0) and np.all(actual[:, :, -1] == 0)
    for name, before in inputs.items():
        np.testing.assert_array_equal(volumes[name], before)


def test_rayleigh_class_energy_rejects_bad_array_shapes() -> None:
    volumes, params = _fixture()
    with pytest.raises(ValueError, match="intensity_field must match"):
        toads_rayleigh_intensity_energy(
            volumes["image_values"], **params, intensity_field=volumes["field"][..., :2]
        )
    with pytest.raises(ValueError, match="class_variance_weights"):
        toads_rayleigh_intensity_energy(
            volumes["image_values"], **{**params, "class_variance_weights": np.ones((3, 3), np.float32)}
        )
