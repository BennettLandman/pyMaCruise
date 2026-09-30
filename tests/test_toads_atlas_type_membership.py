"""Java-referenced 13-class topology dispatch and exact-membership composition."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_atlas_type_membership import (
    map_toads_topologies,
    read_toads_atlas_class_types,
    toads_atlas_type_membership_pass,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_atlas_type_membership"
ATLAS = Path(__file__).parents[1] / "data" / "atlases" / "lesionToads-atlas-2012" / "cruise-atlas-12obj-toads2012.txt"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
EXPECTED_TYPES = ("MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF", *("OBJECT",) * 10)


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _from_bits(values: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32).reshape(shape)


def _parameters() -> dict[str, np.ndarray | float]:
    channels, classes = 2, 13
    centroids = np.empty((channels, classes), dtype=np.float32)
    variances = np.empty_like(centroids)
    for c in range(channels):
        for k in range(classes):
            centroids[c, k] = np.float32(
                np.float32(.19) + np.float32(np.float32(.021) * np.float32(k))
                + np.float32(np.float32(.013) * np.float32(c))
            )
            variances[c, k] = np.float32(
                np.float32(.64) + np.float32(np.float32(.037) * np.float32((k + 2 * c) % 7))
            )
    smooth_prior = np.asarray(
        [np.float32(.55) + np.float32(np.float32(.043) * np.float32(k % 8)) for k in range(classes)],
        dtype=np.float32,
    )
    shape_prior = np.asarray(
        [np.float32(.071) + np.float32(np.float32(.006) * np.float32(k % 6)) for k in range(classes)],
        dtype=np.float32,
    )
    similarity = np.zeros((classes, classes), dtype=np.float32)
    for k in range(classes):
        for m in range(classes):
            if k != m:
                similarity[k, m] = np.float32(
                    np.float32(.035) + np.float32(np.float32(.004) * np.float32((2 * k + 3 * m + 1) % 17))
                )
    exponents = np.ones((channels, classes), dtype=np.float32)
    exponents[:, 1:3] = np.asarray([.78, .91], dtype=np.float32)[:, None]
    return {
        "centroids": centroids,
        "exponent_by_channel_class": exponents,
        "modality_weights": np.asarray([.67, .43], dtype=np.float32),
        "class_variance_weights": variances,
        "image_variances": np.asarray([.58, .37], dtype=np.float32),
        "smoothing_priors": smooth_prior,
        "centroid_similarity": similarity,
        "atlas_shape_priors": shape_prior,
        "smoothing": .039,
        "prior_factor": .62,
        "prior_coefficient": .14,
    }


def _fixture() -> tuple[dict[str, np.ndarray], dict, tuple[dict, ...]]:
    names = (
        "images", "intensity_field", "psf_bias", "initial_memberships",
        "transformed_shapes", "closest_labels", "folded_ribbon", "atlas_template_shapes",
    )
    volumes = {name: nib.load(str(ROOT / f"{name}.nii.gz")) for name in names}
    image = volumes["images"]
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert image.shape == (9, 8, 7, 2)
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (.9, 1.2, 1.45), atol=1e-6, rtol=0)
    for name, volume in volumes.items():
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S"), name
        if name != "atlas_template_shapes":
            np.testing.assert_allclose(volume.affine, image.affine, atol=1e-7, rtol=0)
    atlas = volumes["atlas_template_shapes"]
    assert atlas.shape == (11, 10, 9, 13)
    assert atlas.shape[:3] != image.shape[:3]
    np.testing.assert_allclose(
        atlas.affine @ np.asarray([1, 1, 1, 1], dtype=float),
        image.affine @ np.asarray([0, 0, 0, 1], dtype=float),
        atol=2e-6,
        rtol=0,
    )
    data = {name: np.asarray(volume.dataobj) for name, volume in volumes.items()}
    for name in ("images", "intensity_field", "psf_bias", "initial_memberships", "transformed_shapes", "atlas_template_shapes"):
        data[name] = np.asarray(data[name], dtype=np.float32)
    data["closest_labels"] = np.asarray(data["closest_labels"], dtype=np.int16)
    data["folded_ribbon"] = np.asarray(data["folded_ribbon"], dtype=np.uint8)
    assert data["folded_ribbon"].sum() > 20
    assert not np.array_equal(data["folded_ribbon"], data["folded_ribbon"][::-1])
    return data, _parameters(), tuple(JAVA["cases"])


def _run(data: dict, parameters: dict, case: dict):
    return toads_atlas_type_membership_pass(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        class_types=EXPECTED_TYPES,
        transformed_shapes=data["transformed_shapes"],
        closest_labels=data["closest_labels"],
        intensity_field=data["intensity_field"] if case["field"] else None,
        psf_bias=data["psf_bias"] if case["psf"] else None,
        smoothing=parameters["smoothing"] if case["smoothing"] else 0,
        prior_coefficient=parameters["prior_coefficient"] if case["prior"] else 0,
        **{key: value for key, value in parameters.items() if key not in ("smoothing", "prior_coefficient")},
    )


def test_selected_descriptor_class_order_mapping_and_fallback() -> None:
    assert read_toads_atlas_class_types(ATLAS) == EXPECTED_TYPES
    assert len(EXPECTED_TYPES) == 13
    assert EXPECTED_TYPES.count("MASK") == 1
    assert EXPECTED_TYPES.count("OPTIMIZED_CSF") == 2
    assert EXPECTED_TYPES.count("OBJECT") == 10
    assert map_toads_topologies(("obj", "mask", "optimized_csf", "OBJ", "unknown")) == (
        "OBJECT", "MASK", "OPTIMIZED_CSF", "OBJECT", "OBJECT"
    )
    assert JAVA["class_types"] == list(EXPECTED_TYPES)
    assert JAVA["segmentation_source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert JAVA["atlas_source_sha256"] == "914839a0885cb26ff3ac22daeeb2603ee1960368b788da3ad343a54053f031db"
    assert JAVA["powertable_source_sha256"] == "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
    assert hashlib.sha256(ATLAS.read_bytes()).hexdigest() == "096eb94f9c814f47535e308b8b0d63eaacbbfc646efc1dd06c8ecd8f0892f006"


def test_13_class_membership_composition_matches_selected_java_bits() -> None:
    data, parameters, cases = _fixture()
    assert len(cases) == 16
    assert np.array_equal(parameters["exponent_by_channel_class"][:, 1], parameters["exponent_by_channel_class"][:, 2])
    # The two CSF labels retain distinct model columns even though their
    # category/channel exponents come from one shared PowerTable per channel.
    assert not np.array_equal(parameters["centroids"][:, 1], parameters["centroids"][:, 2])
    assert not np.array_equal(parameters["class_variance_weights"][:, 1], parameters["class_variance_weights"][:, 2])
    labels = data["closest_labels"]
    for k in range(13):
        direct = np.any(labels[..., :2] == k, axis=-1)
        neighbor = np.any(labels[..., 2:] == k, axis=-1) & ~direct
        prior = ~np.any(labels == k, axis=-1)
        assert np.any(direct[1:-1, 1:-1, 1:-1])
        assert np.any(neighbor[1:-1, 1:-1, 1:-1])
        assert np.any(prior[1:-1, 1:-1, 1:-1])

    before = {key: value.copy() for key, value in data.items()}
    for case in cases:
        result = _run(data, parameters, case)
        assert result.class_types == EXPECTED_TYPES
        assert _bits(result.class_data_energy) == case["class_energy_bits"]
        assert _bits(result.memberships) == case["membership_bits"]
        assert _bits(result.denominator) == case["denominator_bits"]
        assert _bits(result.objective) == [case["objective_bits"]]
        assert result.class_data_energy.dtype == np.float32
        assert result.memberships.dtype == np.float32
        assert result.denominator.dtype == np.float32
        assert np.all(np.isfinite(result.memberships))
        np.testing.assert_array_equal(result.memberships[0], data["initial_memberships"][0])
        np.testing.assert_array_equal(result.memberships[-1], data["initial_memberships"][-1])
        np.testing.assert_array_equal(result.memberships[:, 0], data["initial_memberships"][:, 0])
        np.testing.assert_array_equal(result.memberships[:, -1], data["initial_memberships"][:, -1])
        np.testing.assert_array_equal(result.memberships[:, :, 0], data["initial_memberships"][:, :, 0])
        np.testing.assert_array_equal(result.memberships[:, :, -1], data["initial_memberships"][:, :, -1])
        assert np.all(result.denominator[[0, -1], :, :] == 0)
        assert np.all(result.denominator[:, [0, -1], :] == 0)
        assert np.all(result.denominator[:, :, [0, -1]] == 0)
    for key, value in before.items():
        np.testing.assert_array_equal(data[key], value)


@pytest.mark.parametrize(
    "types",
    [
        ("MASK", "OPTIMIZED_CSF", "OBJECT"),
        ("MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF", *("OBJECT",) * 9, "RAYLEIGH"),
        ("MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF", *("OBJECT",) * 9, "OUTLIER"),
    ],
)
def test_selected_13_class_dispatch_rejects_other_orders_or_families(types: tuple[str, ...]) -> None:
    data, parameters, _ = _fixture()
    with pytest.raises(ValueError, match="class_types"):
        toads_atlas_type_membership_pass(
            image_values=data["images"],
            initial_memberships=np.zeros((9, 8, 7, len(types)), dtype=np.float32),
            class_types=types,
            transformed_shapes=np.ones((9, 8, 7, len(types)), dtype=np.float32),
            closest_labels=np.zeros((9, 8, 7, 4), dtype=np.int16),
            smoothing=0,
            prior_coefficient=0,
            **{key: value for key, value in parameters.items() if key not in ("smoothing", "prior_coefficient")},
        )


def test_duplicate_optimized_category_exponents_must_match() -> None:
    data, parameters, cases = _fixture()
    changed = dict(parameters)
    changed["exponent_by_channel_class"] = parameters["exponent_by_channel_class"].copy()
    changed["exponent_by_channel_class"][0, 2] = np.float32(.79)
    with pytest.raises(ValueError, match="share category/channel exponents"):
        _run(data, changed, cases[0])
