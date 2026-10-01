from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_exact_iteration_slice import toads_default_exact_iteration_slice


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_first_thinning"
REFERENCE = json.loads((ROOT / "exact_membership_slice_reference.json").read_text())


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values, dtype="<f4").tobytes()).hexdigest()


def test_default_exact_iteration_slice_matches_java_composition() -> None:
    spec = json.loads((ROOT / "spec.json").read_text())
    parameters = json.loads((ROOT / "membership_stage_spec.json").read_text())
    image_volume = nib.load(str(ROOT / "subject_images.nii.gz"))
    labels_volume = nib.load(str(ROOT / "initial_segmentation.nii.gz"))
    field_volume = nib.load(str(ROOT / "fixed_inhomogeneity_field.nii.gz"))
    shape_volume = nib.load(str(ROOT / "membership_transformed_shapes.nii.gz"))
    relation_volume = nib.load(str(ROOT / "membership_closest_labels.nii.gz"))
    volumes = (image_volume, labels_volume, field_volume, shape_volume, relation_volume)
    for volume in volumes:
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(volume.affine, image_volume.affine)

    images = np.asarray(image_volume.dataobj, dtype=np.float32)
    labels = np.asarray(labels_volume.dataobj, dtype=np.int16)
    field_3d = np.asarray(field_volume.dataobj, dtype=np.float32)
    field = np.broadcast_to(field_3d[..., None], images.shape).copy()
    shapes = np.asarray(shape_volume.dataobj, dtype=np.float32)
    closest = np.asarray(relation_volume.dataobj, dtype=np.int16)
    initial_centroids = np.asarray(spec["centroids"], dtype=np.float32)

    result = toads_default_exact_iteration_slice(
        images,
        initial_centroids,
        labels,
        np.zeros((11, 11), dtype=bool),
        np.asarray(parameters["category_exponents"], dtype=np.float32),
        np.asarray(parameters["modality_weights"], dtype=np.float32),
        np.asarray(parameters["class_variance_weights"], dtype=np.float32),
        np.asarray(parameters["image_variance_weights"], dtype=np.float32),
        shapes,
        closest,
        np.asarray(parameters["smoothing_priors"], dtype=np.float32),
        np.asarray(parameters["atlas_shape_priors"], dtype=np.float32),
        parameters["smoothing"],
        parameters["prior_factor"],
        parameters["prior_coefficient"],
        intensity_field=field,
        crop_border=spec["crop_border"],
        volume_ratio=spec["volume_ratio"],
        thinning_limit=spec["thinning_limit"],
        centroid_smoothness=parameters["centroid_smoothness"],
        prior_scale=parameters["prior_scale"],
    )

    # This fixture's first thinning removes skeleton voxels only; the source
    # categorical labels are the fixed assigned-class input to the statistics.
    np.testing.assert_array_equal(result.first_thinning.thinning.segmentation, labels)
    assert REFERENCE["source"] == (
        "VandyCRUISE MedicAlgorithmToads08 exact default startup -> "
        "computeExactMemberships selected-method composition"
    )
    assert REFERENCE["segmentation_source_sha256"] == (
        "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    )
    assert REFERENCE["preprocess_source_sha256"] == (
        "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7"
    )
    assert REFERENCE["power_table_source_sha256"] == (
        "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
    )
    assert REFERENCE["shape_xyz_channels"] == [36, 36, 36, 2]
    assert REFERENCE["explicit_image_variance_weights"] == [0.61, 0.34]
    assert REFERENCE["smoothing"] == parameters["smoothing"]
    assert REFERENCE["prior_factor"] == parameters["prior_factor"]
    assert REFERENCE["prior_coefficient"] == parameters["prior_coefficient"]
    assert REFERENCE["prior_scale"] == parameters["prior_scale"]
    np.testing.assert_array_equal(_bits(result.refreshed_centroids), REFERENCE["centroids_bits"])
    np.testing.assert_array_equal(
        _bits(result.class_standard_deviations), REFERENCE["standard_deviation_bits"]
    )
    assert _sha256(result.centroid_similarity) == REFERENCE["similarity_sha256"]
    assert _sha256(result.membership_pass.class_data_energy) == REFERENCE["class_energy_sha256"]
    assert _sha256(result.membership_pass.denominator) == REFERENCE["denominator_sha256"]
    assert _sha256(result.membership_pass.memberships) == REFERENCE["membership_bits_sha256"]
    np.testing.assert_array_equal(_bits(result.membership_pass.objective), [REFERENCE["objective_bits"]])
