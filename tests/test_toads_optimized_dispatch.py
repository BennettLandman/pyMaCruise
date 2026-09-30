from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_optimized_dispatch import toads_optimized_membership_dispatch


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_optimized_dispatch"
REFERENCE = ROOT / "java_reference.json"
TOADS_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
POWERTABLE_SHA = "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
TYPES = ("OBJECT", "OPTIMIZED_CSF", "MASK", "OPTIMIZED_GM", "OPTIMIZED_WM")
EXPONENT_BY_TYPE = {"OPTIMIZED_CSF": 0.8, "OPTIMIZED_GM": 1.2, "OPTIMIZED_WM": 1.6}
CENTROIDS = np.asarray(
    [[0.30, 0.35, 0.40, 0.45, 0.50],
     [0.26, 0.31, 0.36, 0.41, 0.46],
     [0.22, 0.27, 0.32, 0.37, 0.42]],
    dtype=np.float32,
)
VARIANCE_WEIGHTS = np.asarray(
    [[0.8, 0.9, 1.0, 1.1, 1.2],
     [1.1, 0.7, 0.9, 1.3, 0.8],
     [0.6, 1.2, 0.75, 0.95, 1.15]],
    dtype=np.float32,
)
MODALITY_WEIGHTS = np.asarray([0.7, 0.5, 0.3], dtype=np.float32)
IMAGE_VARIANCES = np.asarray([0.55, 0.30, 0.15], dtype=np.float32)
SMOOTHING_PRIORS = np.asarray([0.5, 0.7, 0.9, 1.1, 1.3], dtype=np.float32)
ATLAS_PRIORS = np.asarray([0.12, 0.18, 0.22, 0.25, 0.23], dtype=np.float32)
SIMILARITY = np.asarray(
    [[0, .12, .18, .21, .16], [.15, 0, .13, .19, .17],
     [.18, .14, 0, .11, .20], [.16, .18, .13, 0, .12],
     [.19, .15, .17, .14, 0]],
    dtype=np.float32,
)


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _inputs(use_field: bool, use_psf: bool) -> tuple[dict, np.ndarray]:
    image = nib.load(str(ROOT / "images.nii.gz"))
    ribbon_image = nib.load(str(ROOT / "folded_ribbon.nii.gz"))
    field_image = nib.load(str(ROOT / "field.nii.gz"))
    bias_image = nib.load(str(ROOT / "psf_bias.nii.gz"))
    for volume in (image, ribbon_image, field_image, bias_image):
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
        np.testing.assert_allclose(nib.affines.voxel_sizes(volume.affine), [1.05, 1.2, 1.35], rtol=0, atol=1e-6)
        np.testing.assert_allclose(volume.affine, image.affine, rtol=0, atol=1e-7)
    images = np.asarray(image.dataobj, dtype=np.float32)
    ribbon = np.asarray(ribbon_image.dataobj, dtype=np.uint8)
    fields = np.asarray(field_image.dataobj, dtype=np.float32)
    biases = np.asarray(bias_image.dataobj, dtype=np.float32)
    assert images.shape == (*ribbon.shape, 3)
    assert ribbon.sum() > 150
    assert not np.array_equal(ribbon, ribbon[::-1, :, :])

    nx, ny, nz = ribbon.shape
    classes = len(TYPES)
    initial = np.empty((nx, ny, nz, classes), dtype=np.float32)
    shapes = np.empty_like(initial)
    labels = np.full((nx, ny, nz, 4), -1, dtype=np.int8)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for k in range(classes):
                    n_initial = (x + y + z) % 5 if k == 0 else (x + 2 * y + z + 3 * k) % 7
                    initial[x, y, z, k] = np.float32(
                        np.float32(0.12) + np.float32(np.float32(n_initial) * np.float32(0.017))
                    )
                    n_shape = (x + 2 * y + 3 * z + k) % 5
                    shape = np.float32(np.float32(0.2) + np.float32(np.float32(0.5) * ribbon[x, y, z]))
                    shapes[x, y, z, k] = np.float32(shape + np.float32(np.float32(n_shape) * np.float32(0.01)))
                labels[x, y, z, 0] = (x + y + z) % classes
                if (x + 2 * y + z) % 7 == 0:
                    labels[x, y, z, 1] = (labels[x, y, z, 0] + 2) % classes
                if (x + y + 2 * z) % 3 == 0:
                    labels[x, y, z, 2] = (labels[x, y, z, 0] + 1) % classes
                if (x + z) % 4 == 0:
                    labels[x, y, z, 3] = (labels[x, y, z, 0] + 3) % classes

    exponents = np.ones((3, classes), dtype=np.float32)
    for index, kind in enumerate(TYPES):
        if kind in EXPONENT_BY_TYPE:
            exponents[:, index] = np.float32(EXPONENT_BY_TYPE[kind])
    return {
        "image_values": images,
        "initial_memberships": initial,
        "class_types": TYPES,
        "centroids": CENTROIDS,
        "exponent_by_channel_class": exponents,
        "modality_weights": MODALITY_WEIGHTS,
        "class_variance_weights": VARIANCE_WEIGHTS,
        "image_variances": IMAGE_VARIANCES,
        "transformed_shapes": shapes,
        "closest_labels": labels,
        "smoothing_priors": SMOOTHING_PRIORS,
        "centroid_similarity": SIMILARITY,
        "atlas_shape_priors": ATLAS_PRIORS,
        "smoothing": 0.027,
        "prior_factor": 0.64,
        "prior_coefficient": 0.11,
        "intensity_field": fields if use_field else None,
        "psf_bias": biases if use_psf else None,
    }, initial


def test_interleaved_optimized_and_ordinary_types_match_java_and_powertable() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["toads_source_sha256"] == TOADS_SHA
    assert reference["powertable_source_sha256"] == POWERTABLE_SHA
    assert reference["class_types"] == list(TYPES)
    assert len(reference["cases"]) == 4

    for case in reference["cases"]:
        inputs, initial = _inputs(case["field"], case["psf"])
        inputs_before = {
            name: value.copy() for name, value in inputs.items()
            if isinstance(value, np.ndarray)
        }
        result = toads_optimized_membership_dispatch(**inputs)
        assert _bits(result.memberships) == case["updated_bits"]
        assert _bits(result.objective) == [case["objective_bits"]]
        assert result.memberships.dtype == np.float32
        assert result.objective.dtype == np.float32
        # The selected loop updates only interior voxels; all six faces remain
        # byte-for-byte equal to the supplied membership volume.
        np.testing.assert_array_equal(result.memberships[0], initial[0])
        np.testing.assert_array_equal(result.memberships[-1], initial[-1])
        np.testing.assert_array_equal(result.memberships[:, 0], initial[:, 0])
        np.testing.assert_array_equal(result.memberships[:, -1], initial[:, -1])
        np.testing.assert_array_equal(result.memberships[:, :, 0], initial[:, :, 0])
        np.testing.assert_array_equal(result.memberships[:, :, -1], initial[:, :, -1])
        for name, before in inputs_before.items():
            np.testing.assert_array_equal(inputs[name], before)


@pytest.mark.parametrize("kind", ["OUTLIER", "OUTMASK", "LESION", "RAYLEIGH"])
def test_auxiliary_and_unported_types_reject(kind: str) -> None:
    inputs, _ = _inputs(False, False)
    inputs["class_types"] = ("OBJECT", kind, "MASK", "OPTIMIZED_GM", "OPTIMIZED_WM")
    with pytest.raises(ValueError, match="unsupported optimized-dispatch class types"):
        toads_optimized_membership_dispatch(**inputs)


def test_standard_only_lists_reject_and_class_order_is_required() -> None:
    inputs, _ = _inputs(False, False)
    inputs["class_types"] = ("OBJECT", "MASK", "OBJECT", "MASK", "OBJECT")
    with pytest.raises(ValueError, match="at least one optimized tissue"):
        toads_optimized_membership_dispatch(**inputs)
    inputs, _ = _inputs(False, False)
    inputs["class_types"] = ("OBJECT", "MASK", "OPTIMIZED_CSF")
    with pytest.raises(ValueError, match="every membership class"):
        toads_optimized_membership_dispatch(**inputs)


def test_malformed_optimized_exponent_matrix_rejects() -> None:
    inputs, _ = _inputs(False, False)
    inputs["exponent_by_channel_class"] = np.ones((3, 4), dtype=np.float32)
    with pytest.raises(ValueError, match="exponent_by_channel_class must match centroid dimensions"):
        toads_optimized_membership_dispatch(**inputs)
