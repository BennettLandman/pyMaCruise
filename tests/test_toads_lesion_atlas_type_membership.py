"""Selected-Java comparison for the 13-class lesion-bearing TOADS atlas."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_lesion_atlas_type_membership import (
    LESION_TOADS_CLASS_TYPES,
    map_lesion_toads_topologies,
    read_lesion_toads_atlas_class_types,
    toads_lesion_atlas_type_membership_pass,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_lesion_atlas_type_membership"
ATLAS = Path(__file__).parents[1] / "data" / "atlases" / "lesionToads-atlas-2012" / "cruise-atlas-12obj-lesiontoads2012.txt"
JAVA = json.loads((ROOT / "java_reference.json").read_text())


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _parameters() -> dict:
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
        "lesion_centroid": np.asarray([.41, .36], dtype=np.float32),
        "exponent_by_channel_class": exponents,
        "modality_weights": np.asarray([[.67, .22], [.43, .31]], dtype=np.float32),
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
        "transformed_shapes", "closest_labels", "segmentation",
        "initial_lesion_membership", "folded_ribbon", "atlas_template_shapes",
    )
    volumes = {name: nib.load(str(ROOT / f"{name}.nii.gz")) for name in names}
    subject = volumes["images"]
    assert subject.shape == (9, 8, 7, 2)
    assert nib.aff2axcodes(subject.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(subject.affine), (.9, 1.2, 1.45), rtol=0, atol=1e-6)
    for name, volume in volumes.items():
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S"), name
        if name != "atlas_template_shapes":
            np.testing.assert_allclose(volume.affine, subject.affine, rtol=0, atol=1e-7)
    template = volumes["atlas_template_shapes"]
    assert template.shape == (11, 10, 9, 13)
    assert template.shape[:3] != subject.shape[:3]
    np.testing.assert_allclose(
        template.affine @ np.asarray([1, 1, 1, 1], dtype=float),
        subject.affine @ np.asarray([0, 0, 0, 1], dtype=float), rtol=0, atol=2e-6,
    )
    data = {name: np.asarray(volume.dataobj) for name, volume in volumes.items()}
    for name in ("images", "intensity_field", "psf_bias", "initial_memberships", "transformed_shapes", "initial_lesion_membership", "atlas_template_shapes"):
        data[name] = np.asarray(data[name], dtype=np.float32)
    data["closest_labels"] = np.asarray(data["closest_labels"], dtype=np.int16)
    data["segmentation"] = np.asarray(data["segmentation"], dtype=np.int16)
    data["folded_ribbon"] = np.asarray(data["folded_ribbon"], dtype=np.uint8)
    assert data["folded_ribbon"].sum() > 20
    assert not np.array_equal(data["folded_ribbon"], data["folded_ribbon"][::-1])
    return data, _parameters(), tuple(JAVA["cases"])


def _run(data: dict, parameters: dict, case: dict):
    values = dict(parameters)
    values["smoothing"] = values["smoothing"] if case["smoothing"] else 0.0
    values["prior_coefficient"] = values["prior_coefficient"] if case["prior"] else 0.0
    return toads_lesion_atlas_type_membership_pass(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        initial_lesion_membership=data["initial_lesion_membership"],
        segmentation=data["segmentation"],
        class_types=LESION_TOADS_CLASS_TYPES,
        transformed_shapes=data["transformed_shapes"],
        closest_labels=data["closest_labels"],
        intensity_field=data["intensity_field"] if case["field"] else None,
        psf_bias=data["psf_bias"] if case["psf"] else None,
        **values,
    )


def test_selected_lesion_descriptor_order_and_topology_fallback() -> None:
    assert read_lesion_toads_atlas_class_types(ATLAS) == LESION_TOADS_CLASS_TYPES
    assert len(LESION_TOADS_CLASS_TYPES) == 13
    assert LESION_TOADS_CLASS_TYPES[:3] == ("MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF")
    assert LESION_TOADS_CLASS_TYPES[3:10] == ("OBJECT",) * 7
    assert LESION_TOADS_CLASS_TYPES[10:] == ("LESION", "OBJECT", "OBJECT")
    topology = ["mask", "optimized_csf", "optimized_csf", *("obj",) * 7, "lesion", "obj", "obj"]
    topology[5] = "unknown-topology"  # selected constructor's fallback is OBJECT
    assert map_lesion_toads_topologies(topology) == LESION_TOADS_CLASS_TYPES
    assert JAVA["class_types"] == list(LESION_TOADS_CLASS_TYPES)
    assert JAVA["segmentation_source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert JAVA["atlas_source_sha256"] == "914839a0885cb26ff3ac22daeeb2603ee1960368b788da3ad343a54053f031db"
    assert JAVA["powertable_source_sha256"] == "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
    assert hashlib.sha256(ATLAS.read_bytes()).hexdigest() == "e42901c860901e0c562e01e8800eda6ac43c528ea4c7bdf27218ad8daab870d7"


def test_13_class_lesion_membership_composition_matches_java_bits() -> None:
    data, parameters, cases = _fixture()
    assert len(cases) == 16
    assert np.array_equal(parameters["exponent_by_channel_class"][:, 1], parameters["exponent_by_channel_class"][:, 2])
    assert not np.array_equal(parameters["centroids"][:, 1], parameters["centroids"][:, 2])
    assert not np.array_equal(parameters["class_variance_weights"][:, 1], parameters["class_variance_weights"][:, 2])
    lesion_voxels = data["segmentation"] == 10
    assert np.any(lesion_voxels[1:-1, 1:-1, 1:-1])
    assert np.any(~lesion_voxels[1:-1, 1:-1, 1:-1])
    labels = data["closest_labels"]
    for k in range(13):
        direct = np.any(labels[..., :2] == k, axis=-1)
        neighbor = np.any(labels[..., 2:] == k, axis=-1) & ~direct
        prior = ~np.any(labels == k, axis=-1)
        assert np.any(direct[1:-1, 1:-1, 1:-1])
        assert np.any(neighbor[1:-1, 1:-1, 1:-1])
        assert np.any(prior[1:-1, 1:-1, 1:-1])

    before = {key: value.copy() for key, value in data.items()}
    parameters_before = {
        key: value.copy() for key, value in parameters.items()
        if isinstance(value, np.ndarray)
    }
    for case in cases:
        result = _run(data, parameters, case)
        assert result.class_types == LESION_TOADS_CLASS_TYPES
        assert _bits(result.class_data_terms) == case["class_data_bits"]
        assert _bits(result.lesion_data_terms) == case["lesion_data_bits"]
        assert _bits(result.memberships) == case["membership_bits"]
        assert _bits(result.lesion_membership) == case["lesion_map_bits"]
        assert _bits(result.denominator) == case["denominator_bits"]
        assert _bits(result.objective) == [case["objective_bits"]]
        assert result.class_data_terms.dtype == np.float32
        assert result.lesion_data_terms.dtype == np.float32
        assert result.memberships.dtype == np.float32
        assert result.lesion_membership.dtype == np.float32
        for axis in range(3):
            np.testing.assert_array_equal(np.take(result.memberships, 0, axis=axis), np.take(data["initial_memberships"], 0, axis=axis))
            np.testing.assert_array_equal(np.take(result.memberships, -1, axis=axis), np.take(data["initial_memberships"], -1, axis=axis))
        assert np.all(result.denominator[[0, -1], :, :] == 0)
        assert np.all(result.denominator[:, [0, -1], :] == 0)
        assert np.all(result.denominator[:, :, [0, -1]] == 0)
        np.testing.assert_array_equal(result.lesion_data_terms[..., :10], np.zeros_like(result.lesion_data_terms[..., :10]))
        np.testing.assert_array_equal(result.lesion_data_terms[..., 11:], np.zeros_like(result.lesion_data_terms[..., 11:]))
    for key, value in before.items():
        np.testing.assert_array_equal(data[key], value)
    for key, value in parameters_before.items():
        np.testing.assert_array_equal(parameters[key], value)


@pytest.mark.parametrize(
    "types",
    [
        ("MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF", *("OBJECT",) * 10),
        ("MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF", *("OBJECT",) * 6, "LESION", *("OBJECT",) * 3),
    ],
)
def test_other_atlas_orders_reject(types: tuple[str, ...]) -> None:
    data, parameters, _ = _fixture()
    with pytest.raises(ValueError, match="class_types"):
        toads_lesion_atlas_type_membership_pass(
            image_values=data["images"],
            initial_memberships=np.zeros((9, 8, 7, 13), dtype=np.float32),
            initial_lesion_membership=data["initial_lesion_membership"],
            segmentation=data["segmentation"],
            class_types=types,
            transformed_shapes=np.ones((9, 8, 7, 13), dtype=np.float32),
            closest_labels=data["closest_labels"],
            **{key: value for key, value in parameters.items() if key not in ("smoothing", "prior_coefficient")},
            smoothing=0,
            prior_coefficient=0,
        )


def test_category_duplicate_exponents_must_match() -> None:
    data, parameters, cases = _fixture()
    altered = dict(parameters)
    altered["exponent_by_channel_class"] = parameters["exponent_by_channel_class"].copy()
    altered["exponent_by_channel_class"][0, 2] = np.float32(.79)
    with pytest.raises(ValueError, match="share category/channel exponents"):
        _run(data, altered, cases[0])
