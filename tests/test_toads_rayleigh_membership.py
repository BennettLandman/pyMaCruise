from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_rayleigh_membership import toads_rayleigh_membership_pass


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_rayleigh_membership"
SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _from_bits(bits: list[int]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32)


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _fixture() -> tuple[dict, dict, dict]:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    loaded = {name: nib.load(str(ROOT / f"{name}.nii.gz")) for name in (
        "images", "field", "psf_bias", "folded_ribbon", "initial_memberships",
        "transformed_shapes", "closest_labels",
    )}
    image = loaded["images"]
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), [0.9, 1.2, 1.45], rtol=0, atol=1e-6)
    for volume in loaded.values():
        np.testing.assert_allclose(volume.affine, image.affine, rtol=0, atol=1e-7)
    data = {name: np.asarray(volume.dataobj) for name, volume in loaded.items()}
    for name in ("images", "field", "psf_bias", "initial_memberships", "transformed_shapes"):
        data[name] = np.asarray(data[name], dtype=np.float32)
    data["closest_labels"] = np.asarray(data["closest_labels"], dtype=np.int16)
    data["folded_ribbon"] = np.asarray(data["folded_ribbon"], dtype=np.uint8)
    assert data["images"].shape == data["field"].shape == data["psf_bias"].shape == (9, 8, 7, 3)
    assert data["initial_memberships"].shape == data["transformed_shapes"].shape == (9, 8, 7, 3)
    assert data["closest_labels"].shape == (9, 8, 7, 4)
    assert data["folded_ribbon"].sum() > 20

    params = {
        "centroids": _from_bits(reference["centroid_bits"]).reshape(3, 3),
        "modality_weights": _from_bits(reference["modality_bits"]),
        "class_variance_weights": _from_bits(reference["variance_bits"]).reshape(3, 3),
        "image_variances": _from_bits(reference["image_variance_bits"]),
        "smoothing_priors": _from_bits(reference["smoothing_prior_bits"]),
        "atlas_shape_priors": _from_bits(reference["atlas_prior_bits"]),
        "centroid_similarity": _from_bits(reference["similarity_bits"]).reshape(3, 3),
        "smoothing": reference["smoothing"],
        "prior_factor": reference["prior_factor"],
        "prior_coefficient": reference["prior_coefficient"],
    }
    return data, params, reference


def test_rayleigh_membership_composition_matches_selected_java_bits() -> None:
    data, params, reference = _fixture()
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["shape"] == [9, 8, 7, 3]

    labels = data["closest_labels"]
    for class_index in range(3):
        direct = np.any((labels[..., :2] == class_index), axis=-1)
        neighbor = np.any((labels[..., 2:] == class_index), axis=-1) & ~direct
        prior = ~np.any(labels == class_index, axis=-1)
        assert np.any(direct[1:-1, 1:-1, 1:-1])
        assert np.any(neighbor[1:-1, 1:-1, 1:-1])
        assert np.any(prior[1:-1, 1:-1, 1:-1])

    data_before = {name: value.copy() for name, value in data.items()}
    for case in reference["cases"]:
        corrected = (
            data["field"] * data["images"] if case["field"] else data["images"]
        ) - (data["psf_bias"] if case["bias"] else np.float32(0.0))
        interior_values = corrected[1:-1, 1:-1, 1:-1]
        assert np.any(interior_values > 0)
        if case["bias"]:
            assert np.any(interior_values <= 0)
        result = toads_rayleigh_membership_pass(
            image_values=data["images"],
            initial_memberships=data["initial_memberships"],
            transformed_shapes=data["transformed_shapes"],
            closest_labels=data["closest_labels"],
            **params,
            intensity_field=data["field"] if case["field"] else None,
            psf_bias=data["psf_bias"] if case["bias"] else None,
        )
        assert _bits(result.data_energy) == case["energy_bits"]
        assert _bits(result.memberships) == case["membership_bits"]
        assert _bits(result.objective) == [case["objective_bits"]]
        assert result.data_energy.shape == result.memberships.shape == (9, 8, 7, 3)
        assert result.data_energy.dtype == result.memberships.dtype == np.float32
        assert result.objective.dtype == np.float32
        assert np.all(result.memberships[0] == data["initial_memberships"][0])
        assert np.all(result.memberships[-1] == data["initial_memberships"][-1])
        assert np.all(result.memberships[:, 0] == data["initial_memberships"][:, 0])
        assert np.all(result.memberships[:, -1] == data["initial_memberships"][:, -1])
        assert np.all(result.memberships[:, :, 0] == data["initial_memberships"][:, :, 0])
        assert np.all(result.memberships[:, :, -1] == data["initial_memberships"][:, :, -1])
    for name, before in data_before.items():
        np.testing.assert_array_equal(data[name], before)


def test_rayleigh_membership_composition_rejects_nonpositive_scales() -> None:
    data, params, _ = _fixture()
    bad = dict(params)
    bad["centroids"] = params["centroids"].copy()
    bad["centroids"][0, 0] = 0.0
    with pytest.raises(ValueError, match="RAYLEIGH scales must be finite and positive"):
        toads_rayleigh_membership_pass(
            image_values=data["images"],
            initial_memberships=data["initial_memberships"],
            transformed_shapes=data["transformed_shapes"],
            closest_labels=data["closest_labels"],
            **bad,
            intensity_field=data["field"],
            psf_bias=data["psf_bias"],
        )
