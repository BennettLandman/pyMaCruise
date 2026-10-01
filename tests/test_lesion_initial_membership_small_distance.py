from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_initial_memberships import lesion_initial_memberships


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_initial_membership_small_distance"
HARNESS = Path(__file__).parent / "java_reference" / "LesionInitialMembershipSmallDistanceReferenceHarness.java"
SOURCE_SHA = "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08"


def test_initial_membership_zero_and_nonzero_small_distance_threshold_matches_java() -> None:
    metadata = json.loads((ROOT / "inputs.json").read_text())
    reference = json.loads((ROOT / "java_reference.json").read_text())
    image = nib.load(str(ROOT / metadata["nifti"]))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    values = np.asarray(image.dataobj, dtype=np.float32)
    channels_first = np.transpose(values, (3, 0, 1, 2))
    centroids = np.asarray(metadata["centroids_channels_classes"], dtype=np.float32)
    bounds = tuple(tuple(pair) for pair in metadata["crop_bounds_source_coordinates"])
    images_before, centroids_before = channels_first.copy(), centroids.copy()

    result = lesion_initial_memberships(
        channels_first,
        centroids,
        crop_bounds=bounds,
        pad_border=metadata["pad_border"],
    )
    expected = np.asarray(reference["membership_bits"], dtype=np.int32).view(np.float32)
    expected = expected.reshape((*reference["shape"], centroids.shape[1]))
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    assert reference["centroid_bits"] == [int(value.view(np.int32)) for value in centroids.ravel()]
    np.testing.assert_array_equal(result.memberships.view(np.int32), expected.view(np.int32))
    assert result.active_bounds == ((1, 3), (1, 3), (1, 3))
    # This active voxel has exact-zero and nonzero-but-sub-1e-20f squared
    # distances, alongside an above-threshold class distance.
    np.testing.assert_allclose(result.memberships[1, 1, 1], [0.4975, 0.4975, 0.005], atol=2.0e-4)
    for axis_slice in (result.memberships[0], result.memberships[-1]):
        assert np.count_nonzero(axis_slice) == 0
    np.testing.assert_array_equal(channels_first, images_before)
    np.testing.assert_array_equal(centroids, centroids_before)
