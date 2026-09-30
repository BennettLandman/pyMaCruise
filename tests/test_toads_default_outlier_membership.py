from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from pymacruise.toads_default_outlier_membership import toads_default_outlier_membership_pass


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_atlas_outlier"
TOADS_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
POWER_SHA = "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"

_CENTROID = np.asarray(
    [[.28,.31,.34,.37,.40,.43,.46,.49,.52,.55,.58],
     [.24,.27,.30,.33,.36,.39,.42,.45,.48,.51,.54]], dtype=np.float32,
)
_EXPONENT = np.asarray(
    [[1,.78,.78,1.16,1.16,1.16,1.16,1.16,1.55,1.55,1.55],
     [1,.91,.91,1.29,1.29,1.29,1.29,1.29,1.73,1.73,1.73]], dtype=np.float32,
)
_VARIANCE = np.asarray(
    [[.81,.93,1.17,.76,.88,1.02,1.14,.69,.96,1.21,.73],
     [1.19,.72,.87,1.08,.79,1.25,.91,1.13,.68,.98,1.16]], dtype=np.float32,
)
_SMOOTH_PRIOR = np.asarray([.65,.71,.84,.91,1.02,1.13,.76,.88,1.05,1.17,.82], dtype=np.float32)
_ATLAS_PRIOR = np.asarray([.10,.08,.09,.11,.12,.09,.10,.11,.07,.06,.07], dtype=np.float32)
_SIMILARITY = np.asarray([
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
], dtype=np.float32)


def _bits(value: np.ndarray | np.float32) -> list[int]:
    return np.asarray(value, dtype=np.float32).ravel().view(np.int32).tolist()


def _read(name: str) -> np.ndarray:
    return np.asarray(nib.load(str(ROOT / f"{name}.nii.gz")).dataobj)


def _insert_column(values: np.ndarray, column: np.ndarray, index: int) -> np.ndarray:
    return np.concatenate((values[..., :index], column[..., None], values[..., index:]), axis=-1)


def _case_inputs(data: dict[str, np.ndarray], case: dict) -> dict:
    index = case["outlier_index"]
    base_initial = np.asarray(data["initial_memberships"], dtype=np.float32)
    base_shapes = np.asarray(data["transformed_shapes"], dtype=np.float32)
    nx, ny, nz = data["folded_ribbon"].shape
    x, y, z = np.ogrid[:nx, :ny, :nz]
    out_initial = np.asarray(
        np.float32(.16) + ((x + 2 * y + z) % 7).astype(np.float32) * np.float32(.014),
        dtype=np.float32,
    )
    ribbon = np.asarray(data["folded_ribbon"], dtype=np.float32)
    out_shape = np.asarray(
        np.float32(.29) + np.float32(.02) * ribbon
        + ((x + 2 * y + z) % 5).astype(np.float32) * np.float32(.004),
        dtype=np.float32,
    )
    initial = _insert_column(base_initial, out_initial, index)
    shapes = _insert_column(base_shapes, out_shape, index)

    labels = np.asarray(data["closest_labels"], dtype=np.int16).copy()
    for slot in range(4):
        source = labels[..., slot].copy()
        labels[..., slot] = np.where(source >= index, source + 1, source)
    pattern = (x + 2 * y + z) % 3
    labels[..., 0] = np.where(pattern == 0, index, labels[..., 0])
    labels[..., 2] = np.where(pattern == 1, index, labels[..., 2])

    centers = _insert_column(
        _CENTROID,
        np.asarray([.63, .57], dtype=np.float32),
        index,
    )
    exponents = _insert_column(_EXPONENT, np.ones(2, dtype=np.float32), index)
    variances = _insert_column(
        _VARIANCE,
        np.asarray([1.07, .94], dtype=np.float32),
        index,
    )
    smooth = np.insert(_SMOOTH_PRIOR, index, np.float32(.97)).astype(np.float32)
    atlas = np.insert(_ATLAS_PRIOR, index, np.float32(.13)).astype(np.float32)
    similarity = np.zeros((12, 12), dtype=np.float32)
    for k in range(12):
        for m in range(12):
            bk = k if k < index else k - 1
            bm = m if m < index else m - 1
            if k == index or m == index:
                similarity[k, m] = (
                    np.float32(0) if k == m
                    else np.float32(.05) + np.float32((2 * k + 3 * m) % 7) * np.float32(.01)
                )
            else:
                similarity[k, m] = _SIMILARITY[bk, bm]

    return {
        "image_values": np.asarray(data["images"], dtype=np.float32),
        "initial_memberships": initial,
        "initial_outlier_membership": np.asarray(data["initial_outlier_membership"], dtype=np.float32),
        "class_types": (*DEFAULT_ATLAS_CLASS_TYPES[:index], "OUTLIER", *DEFAULT_ATLAS_CLASS_TYPES[index:]),
        "centroids": centers,
        "exponent_by_channel_class": exponents,
        "modality_weights": np.asarray([.73,.46], dtype=np.float32),
        "class_variance_weights": variances,
        "image_variances": np.asarray([.61,.34], dtype=np.float32),
        "outlier_value": 1.3,
        "transformed_shapes": shapes,
        "closest_labels": labels,
        "smoothing_priors": smooth,
        "centroid_similarity": similarity,
        "atlas_shape_priors": atlas,
        "smoothing": .043 if case["smoothing"] else 0.0,
        "prior_factor": .64,
        "prior_coefficient": .16 if case["prior"] else 0.0,
        "intensity_field": data["field"] if case["field"] else None,
        "psf_bias": data["psf_bias"] if case["psf"] else None,
    }


def test_default_optimized_plus_outlier_dispatch_matches_java() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["toads_source_sha256"] == TOADS_SHA
    assert reference["powertable_source_sha256"] == POWER_SHA
    volumes = {
        name: nib.load(str(ROOT / f"{name}.nii.gz"))
        for name in (
            "images", "field", "psf_bias", "folded_ribbon", "initial_memberships",
            "transformed_shapes", "closest_labels", "initial_outlier_membership",
        )
    }
    affine = volumes["images"].affine
    assert nib.aff2axcodes(affine) == ("R", "A", "S")
    for volume in volumes.values():
        np.testing.assert_allclose(volume.affine, affine, rtol=0, atol=1e-7)
    data = {name: np.asarray(image.dataobj) for name, image in volumes.items()}
    data["folded_ribbon"] = np.asarray(data["folded_ribbon"], dtype=np.uint8)
    assert data["images"].shape == (9, 8, 7, 2)
    assert data["initial_memberships"].shape == data["transformed_shapes"].shape == (9, 8, 7, 11)
    assert data["folded_ribbon"].sum() > 20
    assert not np.array_equal(data["folded_ribbon"], data["folded_ribbon"][::-1, :, :])

    cases = reference["cases"]
    case_keys = {
        (c["outlier_index"], c["field"], c["psf"], c["smoothing"], c["prior"])
        for c in cases
    }
    assert len(case_keys) == 18
    assert {0, 6, 11} <= {key[0] for key in case_keys}
    assert {(f, p, s, pr) for _, f, p, s, pr in case_keys if _ == 6} == {
        (f, p, s, pr) for f in (False, True) for p in (False, True)
        for s in (False, True) for pr in (False, True)
    }
    expected_by_key = {
        (c["outlier_index"], c["field"], c["psf"], c["smoothing"], c["prior"]): c
        for c in cases
    }
    immutable_data = {name: array.copy() for name, array in data.items()}

    for key, case in expected_by_key.items():
        inputs = _case_inputs(data, case)
        assert list(inputs["class_types"]) == case["types"]
        outlier_index = inputs["class_types"].index("OUTLIER")
        labels = inputs["closest_labels"]
        direct = np.any(labels[..., :2] == outlier_index, axis=-1)
        neighbor = np.any(labels[..., 2:] == outlier_index, axis=-1) & ~direct
        prior = ~np.any(labels == outlier_index, axis=-1)
        assert direct[1:-1, 1:-1, 1:-1].any()
        assert neighbor[1:-1, 1:-1, 1:-1].any()
        assert prior[1:-1, 1:-1, 1:-1].any()
        before_initial = inputs["initial_memberships"].copy()
        before_outlier = inputs["initial_outlier_membership"].copy()
        result = toads_default_outlier_membership_pass(**inputs)
        assert _bits(result.class_data_energy) == case["class_energy_bits"]
        assert _bits(result.outlier_data_energy) == case["outlier_energy_bits"]
        assert _bits(result.memberships) == case["membership_bits"]
        assert _bits(result.outlier_membership) == case["outlier_map_bits"]
        assert _bits(result.denominator) == case["denominator_bits"]
        assert _bits(result.objective) == [case["objective_bits"]]
        np.testing.assert_array_equal(inputs["initial_memberships"], before_initial)
        np.testing.assert_array_equal(inputs["initial_outlier_membership"], before_outlier)
        np.testing.assert_array_equal(result.memberships[[0, -1]], before_initial[[0, -1]])
        np.testing.assert_array_equal(result.memberships[:, [0, -1]], before_initial[:, [0, -1]])
        np.testing.assert_array_equal(result.memberships[:, :, [0, -1]], before_initial[:, :, [0, -1]])
        np.testing.assert_array_equal(result.outlier_membership[[0, -1]], before_outlier[[0, -1]])
        np.testing.assert_array_equal(result.outlier_membership[:, [0, -1]], before_outlier[:, [0, -1]])
        np.testing.assert_array_equal(result.outlier_membership[:, :, [0, -1]], before_outlier[:, :, [0, -1]])

    for name, values in immutable_data.items():
        np.testing.assert_array_equal(data[name], values)


@pytest.mark.parametrize(
    "types",
    [
        DEFAULT_ATLAS_CLASS_TYPES,
        (*DEFAULT_ATLAS_CLASS_TYPES, "OUTLIER", "OUTLIER"),
        ("OBJECT", "OUTLIER", *DEFAULT_ATLAS_CLASS_TYPES[1:]),
        (*DEFAULT_ATLAS_CLASS_TYPES[:2], "OUTMASK", *DEFAULT_ATLAS_CLASS_TYPES[2:], "OUTLIER"),
    ],
)
def test_unselected_default_outlier_families_reject(types: tuple[str, ...]) -> None:
    data = {name: _read(name) for name in ("images", "initial_memberships", "initial_outlier_membership")}
    if len(types) == data["initial_memberships"].shape[-1]:
        initial = data["initial_memberships"]
    else:
        initial = np.zeros((*data["images"].shape[:3], len(types)), dtype=np.float32)
    with pytest.raises(ValueError, match="selected default 11-class order"):
        toads_default_outlier_membership_pass(
            data["images"], initial, data["initial_outlier_membership"], types,
            np.zeros((2, len(types)), dtype=np.float32),
            np.ones((2, len(types)), dtype=np.float32), np.ones(2, dtype=np.float32),
            np.ones((2, len(types)), dtype=np.float32), np.ones(2, dtype=np.float32), 1.3,
            np.ones((*data["images"].shape[:3], len(types)), dtype=np.float32),
            np.full((*data["images"].shape[:3], 4), -1, dtype=np.int16),
            np.ones(len(types), dtype=np.float32), np.zeros((len(types), len(types)), dtype=np.float32),
            np.ones(len(types), dtype=np.float32), .04, .6, .1,
        )
