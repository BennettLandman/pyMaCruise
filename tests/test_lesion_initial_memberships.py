import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_initial_memberships import lesion_initial_memberships


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_initial_memberships"


def _fixture():
    metadata = json.loads((ROOT / "inputs.json").read_text())
    image = nib.load(str(ROOT / metadata["nifti"]))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    values = np.asarray(image.dataobj, dtype=np.float32)
    channels_first = np.transpose(values, (3, 0, 1, 2))
    centroids = np.asarray(metadata["centroids_channels_classes"], dtype=np.float32)
    bounds = tuple(tuple(pair) for pair in metadata["crop_bounds_source_coordinates"])
    return metadata, image, channels_first, centroids, bounds


def test_lesion_initial_memberships_matches_selected_java_harness_bits():
    metadata, image, data, centroids, bounds = _fixture()
    before_data = data.copy()
    before_centroids = centroids.copy()
    result = lesion_initial_memberships(
        data, centroids, crop_bounds=bounds, pad_border=metadata["pad_border"]
    )
    reference = json.loads((ROOT / "java_reference.json").read_text())

    assert reference["source_sha256"] == (
        "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08"
    )
    assert tuple(reference["shape"]) == data.shape[1:]
    assert reference["centroid_bits"] == [
        int(value.view(np.int32)) for value in centroids.ravel()
    ]
    expected = np.asarray(reference["membership_bits"], dtype=np.int32).view(np.float32)
    expected = expected.reshape(data.shape[1:] + (centroids.shape[1],))
    np.testing.assert_array_equal(result.memberships.view(np.int32), expected.view(np.int32))
    assert result.active_bounds == ((1, 6), (1, 7), (1, 8))
    assert np.count_nonzero(result.memberships[0]) == 0
    assert np.count_nonzero(result.memberships[-1]) == 0
    assert result.memberships[1, 1, 1, 0] == np.float32(2.0)
    assert np.array_equal(data, before_data)
    assert np.array_equal(centroids, before_centroids)
    assert image.shape[-1] == centroids.shape[0]


def test_lesion_initial_memberships_rejects_inconsistent_crop_geometry():
    _, _, data, centroids, bounds = _fixture()
    with pytest.raises(ValueError, match="crop extents"):
        lesion_initial_memberships(
            data, centroids, crop_bounds=((2, 8), *bounds[1:]), pad_border=1
        )
