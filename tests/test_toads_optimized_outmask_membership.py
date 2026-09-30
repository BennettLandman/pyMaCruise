from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_optimized_outmask_membership import (
    toads_optimized_outmask_membership_pass,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_optimized_outmask"
TOADS_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
POWER_SHA = "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
TYPES = ("OPTIMIZED_GM", "OUTMASK", "MASK", "OBJECT", "OPTIMIZED_CSF", "OPTIMIZED_WM")
EXPONENTS = {"OPTIMIZED_CSF": (0.78, 0.91), "OPTIMIZED_GM": (1.16, 1.29), "OPTIMIZED_WM": (1.55, 1.73)}


def _values(bits: list[int]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32)


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _fixture() -> tuple[dict, dict, list[dict]]:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    names = (
        "images", "field", "psf_bias", "folded_ribbon", "initial_memberships",
        "transformed_shapes", "closest_labels", "initial_outlier_membership",
    )
    volumes = {name: nib.load(str(ROOT / f"{name}.nii.gz")) for name in names}
    image = volumes["images"]
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), [0.9, 1.2, 1.45], rtol=0, atol=1e-6)
    for volume in volumes.values():
        np.testing.assert_allclose(volume.affine, image.affine, rtol=0, atol=1e-7)
    data = {name: np.asarray(volume.dataobj) for name, volume in volumes.items()}
    for name in ("images", "field", "psf_bias", "initial_memberships", "transformed_shapes", "initial_outlier_membership"):
        data[name] = np.asarray(data[name], dtype=np.float32)
    data["closest_labels"] = np.asarray(data["closest_labels"], dtype=np.int16)
    data["folded_ribbon"] = np.asarray(data["folded_ribbon"], dtype=np.uint8)
    assert data["images"].shape == (9, 8, 7, 2)
    assert data["initial_memberships"].shape == (9, 8, 7, 6)
    assert data["closest_labels"].shape == (9, 8, 7, 4)
    assert data["folded_ribbon"].sum() > 20
    assert not np.array_equal(data["folded_ribbon"], data["folded_ribbon"][::-1, :, :])
    parameters = {
        "centroids": np.asarray(
            [[.30, .34, .38, .42, .46, .50], [.26, .30, .34, .38, .42, .46]], dtype=np.float32
        ),
        "exponent_by_channel_class": np.asarray(
            [[1.16, 1, 1, 1, .78, 1.55], [1.29, 1, 1, 1, .91, 1.73]], dtype=np.float32
        ),
        "modality_weights": np.asarray([.73, .46], dtype=np.float32),
        "class_variance_weights": np.asarray(
            [[.82, 1.17, .69, .95, 1.23, .76], [1.12, .77, .91, 1.04, .68, 1.31]], dtype=np.float32
        ),
        "image_variances": np.asarray([.61, .34], dtype=np.float32),
        "outlier_value": 1.37,
        "smoothing_priors": np.asarray([.65, .91, 1.13, .72, .87, 1.05], dtype=np.float32),
        "atlas_shape_priors": np.asarray([.15, .21, .12, .19, .16, .17], dtype=np.float32),
        "centroid_similarity": np.asarray(
            [[0,.17,.24,.11,.19,.14],[.13,0,.19,.20,.18,.15],[.22,.16,0,.12,.20,.17],
             [.16,.18,.13,0,.12,.16],[.18,.14,.19,.17,0,.11],[.19,.15,.17,.14,.13,0]],
            dtype=np.float32,
        ),
        "prior_factor": .64,
    }
    return data, parameters, reference["cases"]


def _run(data: dict, parameters: dict, case: dict):
    return toads_optimized_outmask_membership_pass(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        initial_outlier_membership=data["initial_outlier_membership"],
        class_types=TYPES,
        transformed_shapes=data["transformed_shapes"],
        closest_labels=data["closest_labels"],
        smoothing=.043 if case["smoothing"] else 0.0,
        prior_coefficient=.16 if case["prior"] else 0.0,
        intensity_field=data["field"] if case["field"] else None,
        psf_bias=data["psf_bias"] if case["psf"] else None,
        **parameters,
    )


def test_optimized_tissue_outmask_composition_matches_selected_java_bits() -> None:
    data, parameters, cases = _fixture()
    assert len(cases) == 16
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["toads_source_sha256"] == TOADS_SHA
    assert reference["powertable_source_sha256"] == POWER_SHA
    assert reference["class_types"] == list(TYPES)
    labels = data["closest_labels"]
    for k in range(len(TYPES)):
        direct = np.any(labels[..., :2] == k, axis=-1)
        neighbor = np.any(labels[..., 2:] == k, axis=-1) & ~direct
        prior = ~np.any(labels == k, axis=-1)
        assert np.any(direct[1:-1, 1:-1, 1:-1])
        assert np.any(neighbor[1:-1, 1:-1, 1:-1])
        assert np.any(prior[1:-1, 1:-1, 1:-1])

    original = {name: value.copy() for name, value in data.items()}
    result_by_mode = {}
    for case in cases:
        result = _run(data, parameters, case)
        key = (case["field"], case["psf"], case["smoothing"], case["prior"])
        result_by_mode[key] = result
        assert _bits(result.class_data_energy) == case["class_energy_bits"]
        assert _bits(result.outlier_data_energy) == case["outlier_energy_bits"]
        assert _bits(result.memberships) == case["membership_bits"]
        assert _bits(result.outlier_membership) == case["outlier_membership_bits"]
        assert _bits(result.denominator) == case["denominator_bits"]
        assert _bits(result.objective) == [case["objective_bits"]]
        assert result.class_data_energy.dtype == np.float32
        assert result.outlier_data_energy.dtype == np.float32
        assert result.memberships.dtype == np.float32
        assert result.outlier_membership.dtype == np.float32
        assert result.denominator.dtype == np.float32
        np.testing.assert_array_equal(result.memberships[0], data["initial_memberships"][0])
        np.testing.assert_array_equal(result.memberships[-1], data["initial_memberships"][-1])
        np.testing.assert_array_equal(result.memberships[:, 0], data["initial_memberships"][:, 0])
        np.testing.assert_array_equal(result.memberships[:, -1], data["initial_memberships"][:, -1])
        np.testing.assert_array_equal(result.memberships[:, :, 0], data["initial_memberships"][:, :, 0])
        np.testing.assert_array_equal(result.memberships[:, :, -1], data["initial_memberships"][:, :, -1])
        np.testing.assert_array_equal(result.outlier_membership[0], data["initial_outlier_membership"][0])
        np.testing.assert_array_equal(result.outlier_membership[-1], data["initial_outlier_membership"][-1])
        np.testing.assert_array_equal(result.outlier_membership[:, 0], data["initial_outlier_membership"][:, 0])
        np.testing.assert_array_equal(result.outlier_membership[:, -1], data["initial_outlier_membership"][:, -1])
        np.testing.assert_array_equal(result.outlier_membership[:, :, 0], data["initial_outlier_membership"][:, :, 0])
        np.testing.assert_array_equal(result.outlier_membership[:, :, -1], data["initial_outlier_membership"][:, :, -1])
        assert np.all(result.denominator[[0, -1], :, :] == 0)
        assert np.all(result.denominator[:, [0, -1], :] == 0)
        assert np.all(result.denominator[:, :, [0, -1]] == 0)
    for name, before in original.items():
        np.testing.assert_array_equal(data[name], before)

    # OUTMASK data/auxiliary terms ignore field and PSF corrections; optimized
    # categories and OBJECT respond to them where the source does.
    for smoothing, prior in ((False, False), (True, True)):
        a = result_by_mode[(False, False, smoothing, prior)].class_data_energy
        f = result_by_mode[(True, False, smoothing, prior)].class_data_energy
        p = result_by_mode[(False, True, smoothing, prior)].class_data_energy
        np.testing.assert_array_equal(a[..., 1], f[..., 1])
        np.testing.assert_array_equal(a[..., 1], p[..., 1])
        assert not np.array_equal(a[1:-1, 1:-1, 1:-1, 0], f[1:-1, 1:-1, 1:-1, 0])
        assert not np.array_equal(a[1:-1, 1:-1, 1:-1, 0], p[1:-1, 1:-1, 1:-1, 0])
        np.testing.assert_array_equal(
            result_by_mode[(False, False, smoothing, prior)].outlier_data_energy,
            result_by_mode[(True, True, smoothing, prior)].outlier_data_energy,
        )


@pytest.mark.parametrize(
    "types",
    [
        ("OBJECT", "OUTMASK", "MASK"),
        ("OBJECT", "OUTMASK", "OUTMASK", "OPTIMIZED_GM"),
        ("OUTMASK", "OPTIMIZED_GM"),
        ("OBJECT", "MASK", "OUTMASK", "OPTIMIZED_GM", "LESION"),
        ("OBJECT", "MASK", "OUTMASK", "RAYLEIGH"),
        ("OBJECT", "MASK", "OUTMASK", "OPTIMIZED_CSF", "OUTLIER"),
    ],
)
def test_invalid_family_compositions_reject(types: tuple[str, ...]) -> None:
    data, parameters, cases = _fixture()
    case = cases[0]
    with pytest.raises(ValueError, match="class_types must contain"):
        toads_optimized_outmask_membership_pass(
            image_values=data["images"],
            initial_memberships=np.zeros((9, 8, 7, len(types)), dtype=np.float32),
            initial_outlier_membership=data["initial_outlier_membership"],
            class_types=types,
            transformed_shapes=np.ones((9, 8, 7, len(types)), dtype=np.float32),
            closest_labels=np.zeros((9, 8, 7, 4), dtype=np.int16),
            smoothing=0,
            prior_coefficient=0,
            **parameters,
        )
