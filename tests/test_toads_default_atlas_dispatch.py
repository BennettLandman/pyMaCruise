from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_default_atlas_membership import (
    DEFAULT_ATLAS_CLASS_TYPES,
    toads_default_atlas_membership_pass,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_atlas_dispatch"
TOADS_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
POWER_SHA = "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
EXPONENT = {"OPTIMIZED_CSF": (.78, .91), "OPTIMIZED_GM": (1.16, 1.29), "OPTIMIZED_WM": (1.55, 1.73)}


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _values(bits: list[int]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32)


def _fixture() -> tuple[dict, dict, list[dict], dict]:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    names = (
        "images", "field", "psf_bias", "folded_ribbon", "initial_memberships",
        "transformed_shapes", "closest_labels",
    )
    volumes = {name: nib.load(str(ROOT / f"{name}.nii.gz")) for name in names}
    image = volumes["images"]
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), [.9, 1.2, 1.45], rtol=0, atol=1e-6)
    for volume in volumes.values():
        np.testing.assert_allclose(volume.affine, image.affine, rtol=0, atol=1e-7)
    data = {name: np.asarray(volume.dataobj) for name, volume in volumes.items()}
    for name in ("images", "field", "psf_bias", "initial_memberships", "transformed_shapes"):
        data[name] = np.asarray(data[name], dtype=np.float32)
    data["closest_labels"] = np.asarray(data["closest_labels"], dtype=np.int16)
    data["folded_ribbon"] = np.asarray(data["folded_ribbon"], dtype=np.uint8)
    assert data["images"].shape == (9, 8, 7, 2)
    assert data["initial_memberships"].shape == data["transformed_shapes"].shape == (9, 8, 7, 11)
    assert data["folded_ribbon"].sum() > 20
    assert not np.array_equal(data["folded_ribbon"], data["folded_ribbon"][::-1, :, :])

    types = DEFAULT_ATLAS_CLASS_TYPES
    parameters = {
        "class_types": types,
        "centroids": np.asarray(
            [[.28,.31,.34,.37,.40,.43,.46,.49,.52,.55,.58], [.24,.27,.30,.33,.36,.39,.42,.45,.48,.51,.54]],
            dtype=np.float32,
        ),
        "exponent_by_channel_class": np.asarray(
            [[1,.78,.78,1.16,1.16,1.16,1.16,1.16,1.55,1.55,1.55],
             [1,.91,.91,1.29,1.29,1.29,1.29,1.29,1.73,1.73,1.73]], dtype=np.float32,
        ),
        "modality_weights": np.asarray([.73,.46], dtype=np.float32),
        "class_variance_weights": np.asarray(
            [[.81,.93,1.17,.76,.88,1.02,1.14,.69,.96,1.21,.73],
             [1.19,.72,.87,1.08,.79,1.25,.91,1.13,.68,.98,1.16]], dtype=np.float32,
        ),
        "image_variances": np.asarray([.61,.34], dtype=np.float32),
        "smoothing_priors": np.asarray([.65,.71,.84,.91,1.02,1.13,.76,.88,1.05,1.17,.82], dtype=np.float32),
        "atlas_shape_priors": np.asarray([.10,.08,.09,.11,.12,.09,.10,.11,.07,.06,.07], dtype=np.float32),
        "centroid_similarity": np.asarray([
            [0,.11,.14,.18,.13,.16,.12,.15,.09,.10,.08],
            [.12,0,.13,.17,.11,.16,.14,.10,.09,.12,.08],
            [.15,.12,0,.13,.18,.12,.16,.11,.10,.07,.09],
            [.14,.16,.12,0,.13,.18,.11,.14,.08,.10,.12],
            [.11,.13,.17,.12,0,.11,.19,.14,.10,.09,.07],
            [.16,.12,.14,.17,.11,0,.12,.18,.09,.08,.10],
            [.13,.15,.11,.14,.18,.12,0,.13,.08,.09,.11],
            [.12,.10,.14,.16,.13,.18,.11,0,.07,.09,.08],
            [.10,.09,.12,.08,.11,.09,.10,.07,0,.13,.16],
            [.11,.12,.08,.10,.09,.07,.11,.08,.14,0,.12],
            [.08,.10,.09,.12,.07,.11,.10,.09,.17,.13,0],
        ], dtype=np.float32),
        "prior_factor": .64,
    }
    assert reference["class_types"] == list(types)
    return data, parameters, reference["cases"], reference


def _run(data: dict, parameters: dict, case: dict):
    return toads_default_atlas_membership_pass(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        transformed_shapes=data["transformed_shapes"],
        closest_labels=data["closest_labels"],
        smoothing=.043 if case["smoothing"] else 0.0,
        prior_coefficient=.16 if case["prior"] else 0.0,
        intensity_field=data["field"] if case["field"] else None,
        psf_bias=data["psf_bias"] if case["psf"] else None,
        **parameters,
    )


def test_default_wrapper_ordered_11_class_dispatch_matches_java() -> None:
    data, parameters, cases, reference = _fixture()
    assert reference["toads_source_sha256"] == TOADS_SHA
    assert reference["powertable_source_sha256"] == POWER_SHA
    assert DEFAULT_ATLAS_CLASS_TYPES == (
        "MASK", "OPTIMIZED_CSF", "OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_GM",
        "OPTIMIZED_GM", "OPTIMIZED_GM", "OPTIMIZED_GM", "OPTIMIZED_WM",
        "OPTIMIZED_WM", "OPTIMIZED_WM",
    )
    assert len(cases) == 16

    labels = data["closest_labels"]
    for k in range(11):
        direct = np.any(labels[..., :2] == k, axis=-1)
        neighbor = np.any(labels[..., 2:] == k, axis=-1) & ~direct
        prior = ~np.any(labels == k, axis=-1)
        for route in (direct, neighbor, prior):
            assert route[1:-1, 1:-1, 1:-1].any(), f"class {k} misses a relation route"

    source_inputs = {name: value.copy() for name, value in data.items()}
    result_by_mode = {}
    for case in cases:
        result = _run(data, parameters, case)
        key = (case["field"], case["psf"], case["smoothing"], case["prior"])
        result_by_mode[key] = result
        assert _bits(result.class_data_energy) == case["class_energy_bits"]
        assert _bits(result.memberships) == case["membership_bits"]
        assert _bits(result.denominator) == case["denominator_bits"]
        assert _bits(result.objective) == [case["objective_bits"]]
        assert result.class_data_energy.dtype == np.float32
        assert result.memberships.dtype == np.float32
        assert result.denominator.dtype == np.float32
        np.testing.assert_array_equal(result.memberships[0], data["initial_memberships"][0])
        np.testing.assert_array_equal(result.memberships[-1], data["initial_memberships"][-1])
        np.testing.assert_array_equal(result.memberships[:, 0], data["initial_memberships"][:, 0])
        np.testing.assert_array_equal(result.memberships[:, -1], data["initial_memberships"][:, -1])
        np.testing.assert_array_equal(result.memberships[:, :, 0], data["initial_memberships"][:, :, 0])
        np.testing.assert_array_equal(result.memberships[:, :, -1], data["initial_memberships"][:, :, -1])
        assert not np.any(result.denominator[0]) and not np.any(result.denominator[-1])

    for name, original in source_inputs.items():
        np.testing.assert_array_equal(data[name], original)

    # Each duplicate optimized category uses its shared category exponent, but
    # retains its own centroid/variance column; field and PSF change those terms.
    for category in ("OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM"):
        indexes = [i for i, kind in enumerate(DEFAULT_ATLAS_CLASS_TYPES) if kind == category]
        assert len({tuple(parameters["centroids"][:, i]) for i in indexes}) == len(indexes)
        for i in indexes:
            assert tuple(parameters["exponent_by_channel_class"][:, i]) == EXPONENT[category]
    for kind in (False, True):
        ordinary = result_by_mode[(False, False, kind, True)].class_data_energy
        with_field = result_by_mode[(True, False, kind, True)].class_data_energy
        with_psf = result_by_mode[(False, True, kind, True)].class_data_energy
        np.testing.assert_array_equal(ordinary[..., 0], with_field[..., 0])
        np.testing.assert_array_equal(ordinary[..., 0], with_psf[..., 0])
        assert not np.array_equal(ordinary[1:-1,1:-1,1:-1,1:], with_field[1:-1,1:-1,1:-1,1:])
        assert not np.array_equal(ordinary[1:-1,1:-1,1:-1,1:], with_psf[1:-1,1:-1,1:-1,1:])

    # Without a positive prior coefficient the selected Java local shape
    # vector remains zero, suppressing classes that are PRIOR-only.
    no_prior = result_by_mode[(False, False, False, False)]
    for k in range(11):
        prior_only = ~np.any(labels == k, axis=-1)
        prior_only[[0, -1], :, :] = False
        prior_only[:, [0, -1], :] = False
        prior_only[:, :, [0, -1]] = False
        assert np.all(no_prior.class_data_energy[..., k][prior_only] == 0)
        assert np.all(no_prior.memberships[..., k][prior_only] == 0)


@pytest.mark.parametrize(
    "types",
    [
        ("MASK",) + ("OPTIMIZED_CSF",) * 2 + ("OPTIMIZED_GM",) * 4 + ("OPTIMIZED_WM",) * 4,
        ("MASK",) + ("OPTIMIZED_CSF",) * 2 + ("OPTIMIZED_GM",) * 5 + ("OPTIMIZED_WM",) * 2 + ("OUTLIER",),
        ("MASK",) + ("OPTIMIZED_CSF",) * 2 + ("OPTIMIZED_GM",) * 5 + ("OPTIMIZED_WM",) * 2 + ("LESION",),
    ],
)
def test_alternate_or_special_type_lists_reject(types: tuple[str, ...]) -> None:
    data, parameters, cases, _ = _fixture()
    with pytest.raises(ValueError, match="selected default 11-class atlas order"):
        toads_default_atlas_membership_pass(
            image_values=data["images"],
            initial_memberships=np.zeros((9, 8, 7, 11), dtype=np.float32),
            transformed_shapes=np.ones((9, 8, 7, 11), dtype=np.float32),
            closest_labels=np.zeros((9, 8, 7, 4), dtype=np.int16),
            smoothing=0,
            prior_coefficient=0,
            class_types=types,
            **{key: value for key, value in parameters.items() if key != "class_types"},
        )


def test_optimized_exponents_must_follow_category_cache() -> None:
    data, parameters, cases, _ = _fixture()
    changed = dict(parameters)
    changed["exponent_by_channel_class"] = parameters["exponent_by_channel_class"].copy()
    changed["exponent_by_channel_class"][0, 2] = np.float32(.79)
    with pytest.raises(ValueError, match="share category/channel exponents"):
        toads_default_atlas_membership_pass(
            image_values=data["images"],
            initial_memberships=data["initial_memberships"],
            transformed_shapes=data["transformed_shapes"],
            closest_labels=data["closest_labels"],
            smoothing=0,
            prior_coefficient=0,
            **changed,
        )
