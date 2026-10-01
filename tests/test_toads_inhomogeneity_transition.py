from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_inhomogeneity_state import ToadsInhomogeneityFieldState
from pymacruise.toads_inhomogeneity_transition import (
    toads_corrected_channel_samples,
    toads_next_inhomogeneity_degree,
    update_toads_approximate_field_transition,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_inhomogeneity_transition"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsInhomogeneityTransitionReferenceHarness.java"
HARNESS_SHA256 = "6d33a7d66b8e0560ea1f525182304499e0f9f984f14303ee39da1a60ed8f0794"
SOURCE_HASHES = {
    "algorithm_sha256": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
    "segmentation_sha256": "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884",
    "correction_sha256": "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8",
    "jama_wrapper_sha256": "efd52d4a8f86984c1dc52d13e01175e07120ccc1e9867283846a0b17f1025526",
    "jama_matrix_sha256": "76ab1241a91518c3b71db14d19a2a3d6c35376f65fdc9804fe26fa79e964ff9c",
    "jama_lu_sha256": "9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6",
    "power_table_sha256": "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678",
}
COEFFICIENT_MAX_ABS = {2: 3.3306690738754696e-13, 3: 1.4509726753431096e-11}


def _load(name: str) -> np.ndarray:
    volume = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
    np.testing.assert_allclose(
        nib.affines.voxel_sizes(volume.affine), (.85, 1.15, 1.55), rtol=0, atol=1e-6
    )
    return np.asarray(volume.dataobj).copy()


def _from_double_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.int64).view(np.float64)


def _field_hash(field: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(field, dtype="<f4").tobytes(order="F")).hexdigest()


def _check_fit(actual, reference: dict, degree: int, field: np.ndarray) -> None:
    assert int(actual.correction_mask.sum()) == reference["mask_count"] == 2802
    assert reference["sample_count"] == 94
    np.testing.assert_array_equal(actual.right_hand_side, _from_double_bits(reference["rhs_bits"]))
    np.testing.assert_array_equal(
        actual.normal_matrix.ravel(), _from_double_bits(reference["normal_matrix_bits"])
    )
    expected_coefficients = _from_double_bits(reference["coefficient_bits"])
    assert float(np.max(np.abs(actual.coefficients - expected_coefficients))) <= COEFFICIENT_MAX_ABS[degree]
    assert _field_hash(field) == reference["field_sha256"]


def test_selected_ratio_gates_and_degree_cap_match_java() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for key, expected in SOURCE_HASHES.items():
        assert JAVA[key] == expected
    assert JAVA["degree_path"] == [2, 3]
    assert JAVA["boundary_degrees"] == [3, 2, 2, 3]
    low = np.float32(JAVA["low_ratio"])
    below = np.nextafter(low, np.float32(-np.inf), dtype=np.float32)
    above = np.nextafter(low, np.float32(np.inf), dtype=np.float32)
    assert toads_next_inhomogeneity_degree(1, 3, below, low) == 2
    assert toads_next_inhomogeneity_degree(2, 3, low, low) == 2
    assert toads_next_inhomogeneity_degree(2, 3, above, low) == 2
    assert toads_next_inhomogeneity_degree(3, 3, below, low) == 3


def test_two_channel_degree2_to3_refit_reuses_attached_prior_field() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for key, expected in SOURCE_HASHES.items():
        assert JAVA[key] == expected
    images = _load("images").astype(np.float32)
    memberships = _load("memberships").astype(np.float32)
    best = _load("best")
    segmentation = _load("segmentation")
    initial_fields = _load("attached_fields").astype(np.float32)
    centers = np.loadtxt(ROOT / "centroids.txt", dtype=np.float32)
    originals = tuple(x.copy() for x in (images, memberships, best, segmentation))

    field_refs = [initial_fields[..., channel].copy() for channel in range(2)]
    state = ToadsInhomogeneityFieldState(channel_count=2)
    for channel, field in enumerate(field_refs):
        state.add_inhomogeneity_correction(field, channel)
    attached = tuple(state.field_by_channel)
    assert all(attached[channel] is field_refs[channel] for channel in range(2))
    assert JAVA["attached_reference_stable"] is True

    first = update_toads_approximate_field_transition(
        images, memberships, best, segmentation, centers, attached,
        degree=1, polynomial_degree=3, previous_energy_ratio=0.4,
        low_energy_ratio=0.5,
    )
    assert first.degree == 2 and len(first.fits) == 2
    assert all(first.fields[c] is attached[c] for c in range(2))
    assert all(first.fits[c].field is attached[c] for c in range(2))
    for channel in range(2):
        _check_fit(
            first.fits[channel], JAVA["degree2_channel_fits"][channel], 2, attached[channel]
        )
    degree2_fields = tuple(field.copy() for field in attached)

    second = update_toads_approximate_field_transition(
        images, memberships, best, segmentation, centers, attached,
        degree=first.degree, polynomial_degree=3, previous_energy_ratio=0.4,
        low_energy_ratio=0.5,
    )
    assert second.degree == 3 and len(second.fits) == 2
    assert all(second.fields[c] is attached[c] for c in range(2))
    assert all(second.fits[c].field is attached[c] for c in range(2))
    for channel in range(2):
        _check_fit(
            second.fits[channel], JAVA["degree3_channel_fits"][channel], 3, attached[channel]
        )
        # The initial attached field is one everywhere; each source refit only
        # writes the non-background mask and leaves its exterior untouched.
        np.testing.assert_array_equal(degree2_fields[channel][segmentation == 0], 1.0)
        np.testing.assert_array_equal(attached[channel][segmentation == 0], 1.0)

    assert [
        _field_hash(attached[channel]) for channel in range(2)
    ] == JAVA["final_field_sha256"]
    coordinates = np.asarray(JAVA["corrected_sample_coordinates"], dtype=np.int64)
    corrected = toads_corrected_channel_samples(images, np.stack(attached, axis=-1), coordinates)
    expected_bits = np.asarray(JAVA["corrected_sample_bits"], dtype=np.int32).view(np.float32)
    np.testing.assert_array_equal(corrected.view(np.int32), expected_bits.view(np.int32))
    for current, original in zip((images, memberships, best, segmentation), originals, strict=True):
        np.testing.assert_array_equal(current, original)
