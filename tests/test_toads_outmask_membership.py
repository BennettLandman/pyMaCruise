from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_outmask_membership import toads_outmask_membership_pass


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_outmask_membership"
SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _float(bits: int) -> np.float32:
    return np.asarray([bits], dtype=np.int32).view(np.float32)[0]


def _values(bits: list[int]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32)


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _fixture() -> tuple[dict, dict, dict]:
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
    assert data["images"].shape == data["field"].shape == data["psf_bias"].shape == (9, 8, 7, 2)
    assert data["initial_memberships"].shape == data["transformed_shapes"].shape == (9, 8, 7, 3)
    assert data["closest_labels"].shape == (9, 8, 7, 4)
    assert data["initial_outlier_membership"].shape == (9, 8, 7)
    assert data["folded_ribbon"].sum() > 20
    params = {
        "centroids": _values(reference["centroid_bits"]).reshape(2, 3),
        "modality_weights": _values(reference["modality_bits"]),
        "class_variance_weights": _values(reference["variance_bits"]).reshape(2, 3),
        "image_variances": _values(reference["image_variance_bits"]),
        "smoothing_priors": _values(reference["smoothing_prior_bits"]),
        "atlas_shape_priors": _values(reference["atlas_prior_bits"]),
        "centroid_similarity": _values(reference["similarity_bits"]).reshape(3, 3),
        "outlier_value": float(_float(reference["outlier_value_bits"])),
        "smoothing": reference["smoothing"],
        "prior_factor": reference["prior_factor"],
        "prior_coefficient": reference["prior_coefficient"],
    }
    return data, params, reference


def test_outmask_membership_pass_matches_java_including_aux_smoothing_quirk() -> None:
    data, params, reference = _fixture()
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["class_types"] == ["MASK", "OUTMASK", "OBJECT"]
    labels = data["closest_labels"]
    for k in range(3):
        direct = np.any(labels[..., :2] == k, axis=-1)
        neighbor = np.any(labels[..., 2:] == k, axis=-1) & ~direct
        prior = ~np.any(labels == k, axis=-1)
        assert np.any(direct[1:-1, 1:-1, 1:-1])
        assert np.any(neighbor[1:-1, 1:-1, 1:-1])
        assert np.any(prior[1:-1, 1:-1, 1:-1])

    original = {name: value.copy() for name, value in data.items()}
    for case in reference["cases"]:
        result = toads_outmask_membership_pass(
            image_values=data["images"],
            initial_memberships=data["initial_memberships"],
            initial_outlier_membership=data["initial_outlier_membership"],
            class_types=reference["class_types"],
            transformed_shapes=data["transformed_shapes"],
            closest_labels=data["closest_labels"],
            **params,
            intensity_field=data["field"] if case["field"] else None,
            psf_bias=data["psf_bias"] if case["bias"] else None,
        )
        assert _bits(result.class_data_energy) == case["class_energy_bits"]
        assert _bits(result.outlier_data_energy) == case["outlier_energy_bits"]
        assert _bits(result.memberships) == case["membership_bits"]
        assert _bits(result.outlier_membership) == case["outlier_membership_bits"]
        assert _bits(result.objective) == [case["objective_bits"]]
        assert result.class_data_energy.dtype == result.outlier_data_energy.dtype == np.float32
        assert result.memberships.dtype == result.outlier_membership.dtype == np.float32
        assert result.objective.dtype == np.float32
        assert np.all(result.class_data_energy[0] == 0)
        assert np.all(result.outlier_data_energy[0] == 0)
        assert np.array_equal(result.memberships[0], data["initial_memberships"][0])
        assert np.array_equal(result.outlier_membership[0], data["initial_outlier_membership"][0])
    for name, before in original.items():
        np.testing.assert_array_equal(data[name], before)

    # OUTMASK/MASK and auxiliary data terms ignore field/PSF inputs, while
    # the OBJECT class's energy changes with the selected correction.
    cases = reference["cases"]
    for left, right in ((cases[0], cases[1]), (cases[2], cases[3])):
        left_energy = _values(left["class_energy_bits"]).reshape(9, 8, 7, 3)
        right_energy = _values(right["class_energy_bits"]).reshape(9, 8, 7, 3)
        np.testing.assert_array_equal(left_energy[..., 0], right_energy[..., 0])
        np.testing.assert_array_equal(left_energy[..., 1], right_energy[..., 1])
        assert not np.array_equal(
            left_energy[1:-1, 1:-1, 1:-1, 2], right_energy[1:-1, 1:-1, 1:-1, 2]
        )
    assert cases[0]["outlier_energy_bits"] == cases[1]["outlier_energy_bits"]
    assert cases[2]["outlier_energy_bits"] == cases[3]["outlier_energy_bits"]
    bias_off = _values(cases[0]["class_energy_bits"]).reshape(9, 8, 7, 3)
    field_on = _values(cases[2]["class_energy_bits"]).reshape(9, 8, 7, 3)
    assert not np.array_equal(bias_off[1:-1, 1:-1, 1:-1, 2], field_on[1:-1, 1:-1, 1:-1, 2])


@pytest.mark.parametrize(
    ("types", "match"),
    [(["OBJECT", "MASK"], "exactly one OUTMASK"), (["OUTMASK"], "OBJECT or MASK"),
     (["OBJECT", "OUTMASK", "OUTMASK"], "exactly one OUTMASK"),
     (["OBJECT", "LESION", "OUTMASK"], "OBJECT or MASK")],
)
def test_outmask_class_type_contract(types: list[str], match: str) -> None:
    data, params, _ = _fixture()
    with pytest.raises(ValueError, match=match):
        toads_outmask_membership_pass(
            image_values=data["images"],
            initial_memberships=data["initial_memberships"],
            initial_outlier_membership=data["initial_outlier_membership"],
            class_types=types,
            transformed_shapes=data["transformed_shapes"],
            closest_labels=data["closest_labels"],
            **params,
        )


def test_outmask_auxiliary_input_contract() -> None:
    data, params, reference = _fixture()
    common = {
        "image_values": data["images"],
        "initial_memberships": data["initial_memberships"],
        "class_types": reference["class_types"],
        "transformed_shapes": data["transformed_shapes"],
        "closest_labels": data["closest_labels"],
        **params,
    }
    with pytest.raises(TypeError, match="initial_outlier_membership"):
        toads_outmask_membership_pass(**common)
    with pytest.raises(ValueError, match="initial_outlier_membership must match"):
        toads_outmask_membership_pass(
            **common,
            initial_outlier_membership=data["initial_outlier_membership"][..., None],
        )
    with pytest.raises(ValueError, match="outlier_value must be finite"):
        toads_outmask_membership_pass(
            **{**common, "outlier_value": np.inf},
            initial_outlier_membership=data["initial_outlier_membership"],
        )
