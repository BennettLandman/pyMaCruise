from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_exact_centroids import toads_default_exact_centroid_update
from pymacruise.toads_default_exact_variances import toads_default_exact_class_variances
from pymacruise.toads_first_thinning import toads_first_exact_thinning_pass


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_first_thinning"
CENTROID_REFERENCE = json.loads((ROOT / "centroid_reference.json").read_text())
REFERENCE = json.loads((ROOT / "variance_reference.json").read_text())


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def test_default_exact_assigned_class_variances_match_java_after_centroids() -> None:
    spec = json.loads((ROOT / "spec.json").read_text())
    images_nifti = nib.load(str(ROOT / "subject_images.nii.gz"))
    labels_nifti = nib.load(str(ROOT / "initial_segmentation.nii.gz"))
    field_nifti = nib.load(str(ROOT / "fixed_inhomogeneity_field.nii.gz"))
    for volume in (images_nifti, labels_nifti, field_nifti):
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(volume.affine, images_nifti.affine)
    images = np.asarray(images_nifti.dataobj, dtype=np.float32)[..., :1]
    labels = np.asarray(labels_nifti.dataobj, dtype=np.int16)
    field = np.asarray(field_nifti.dataobj, dtype=np.float32)[..., None]
    initial_centroids = np.asarray(spec["centroids"], dtype=np.float32)[:1]

    initialized = toads_first_exact_thinning_pass(
        np.moveaxis(images, -1, 0),
        initial_centroids,
        labels,
        np.zeros((11, 11), dtype=bool),
        np.zeros((*labels.shape, 11), dtype=np.float32),
        crop_border=spec["crop_border"],
        volume_ratio=spec["volume_ratio"],
        limit=spec["thinning_limit"],
    )
    np.testing.assert_array_equal(initialized.thinning.segmentation, labels)
    refreshed_centroids = toads_default_exact_centroid_update(
        images,
        initialized.memberships,
        initial_centroids,
        np.asarray([[0.8, 1.0, 1.0]], dtype=np.float32),
        intensity_field=field,
        smoothness=0.1,
    )
    np.testing.assert_array_equal(
        _bits(refreshed_centroids), CENTROID_REFERENCE["updated_centroids_bits"]
    )
    actual = toads_default_exact_class_variances(
        images,
        initialized.memberships,
        refreshed_centroids,
        initialized.thinning.segmentation,
        intensity_field=field,
    )

    assert REFERENCE["source"] == (
        "VandyCRUISE MultipleToadSegmentation.computeExactVariances "
        "default MASK+optimized 10obj dispatch"
    )
    assert REFERENCE["segmentation_source_sha256"] == (
        "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    )
    assert REFERENCE["preprocess_source_sha256"] == (
        "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7"
    )
    assert REFERENCE["atlas_descriptor_sha256"] == (
        "cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939"
    )
    assert REFERENCE["shape_xyz_channels"] == [36, 36, 36, 1]
    assert REFERENCE["classes"] == 11
    np.testing.assert_array_equal(_bits(refreshed_centroids), REFERENCE["centroids_bits"])
    np.testing.assert_array_equal(_bits(actual), REFERENCE["standard_deviation_bits"])
