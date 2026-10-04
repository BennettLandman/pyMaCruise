from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_registered_initial_membership import (
    run_lesion_toads_registered_initial_membership,
)
from pymacruise.lesion_brainstem import lesion_toads_brainstem_distance
from pymacruise.lesion_gm_distance import lesion_toads_gm_distance
from pymacruise.lesion_initial_thinning_composition import run_lesion_initial_thinning_pass
from pymacruise.lesion_relations import lesion_compute_relations
from pymacruise.lesion_statistics import (
    compute_lesion_toads_channel_weights,
    compute_lesion_toads_variances,
)
from pymacruise.lesion_standard13_voxel import STANDARD13_NAMES, STANDARD13_TYPES
from pymacruise.lesion_standard13_volume import run_lesion_toads_standard13_volume
from pymacruise.lesion_ventricle import lesion_toads_ventricle_distance


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "lesion_registered_initial_membership"


def _sha256(array: np.ndarray, dtype: str) -> str:
    return hashlib.sha256(array.astype(dtype, copy=False).tobytes()).hexdigest()


def test_registered_initial_membership_and_first_thinning_on_synthetic_ras_grid() -> None:
    spec = json.loads((FIXTURE / "spec.json").read_text())
    image_nii = nib.load(FIXTURE / "images_ras.nii.gz")
    labels_nii = nib.load(FIXTURE / "constructor_labels_ras.nii.gz")
    shapes_nii = nib.load(FIXTURE / "transformed_shapes_ras.nii.gz")
    assert aff2axcodes(image_nii.affine) == aff2axcodes(labels_nii.affine) == aff2axcodes(shapes_nii.affine) == (
        "R",
        "A",
        "S",
    )
    np.testing.assert_array_equal(image_nii.affine, np.asarray(spec["affine_ras_mm"]))
    images = np.asanyarray(image_nii.dataobj).astype(np.float32)
    labels = np.asanyarray(labels_nii.dataobj).astype(np.int16)
    shapes = np.asanyarray(shapes_nii.dataobj).astype(np.float32)
    assert _sha256(images, "<f4") == spec["image_sha256_float32le"]
    assert _sha256(labels, "<i2") == spec["labels_sha256_int16le"]
    assert _sha256(shapes, "<f4") == spec["shapes_sha256_float32le"]

    # These explicit synthetic centroids model the already established
    # pre-registration and post-registration centroid inputs at this seam.
    class_centroids = np.stack(
        (np.linspace(0.2, 0.7, 13), np.linspace(0.25, 0.8, 13))
    ).astype(np.float32)
    result = run_lesion_toads_registered_initial_membership(
        images,
        labels,
        class_centroids,
        np.asarray((0.57, 0.63), dtype=np.float32),
        np.asarray((0.34, 0.4), dtype=np.float32),
        shapes,
        labels.copy(),
        np.zeros((*labels.shape, 13), dtype=np.bool_),
        np.zeros(labels.shape, dtype=np.float32),
        iterations=1,
        first_thinning_limit=1000.0,
    )

    assert result.stage_order[:3] == (
        "constructorState:alreadyOnRegisteredRecropGrid",
        "setGMflag:beforeInitialMemberships",
        "computeRelations",
    )
    assert "computeGMDistance:GMflag=true" in result.stage_order
    assert "computeMemberships:GMflag=true" in result.stage_order
    assert result.excluded_wrapper_stages == (
        "atlas.updateShapeRegistration:transformAndTransformedShapesAreSupplied",
        "InhomogeneityCorrection2:disabledInSelectedMaCRUISEplusProfile",
    )
    assert result.initial_class_memberships.shape == (*labels.shape, 13)
    assert np.count_nonzero(result.initial_class_memberships) == 0
    assert np.all(result.initial_outlier_membership == np.float32(0.5))
    assert np.all(result.initial_lesion_relation == np.float32(1.0))
    # The initialized lesion map is empty; retain and explicitly exercise the
    # source-faithful zero-deviation branch rather than silently regularizing it.
    np.testing.assert_array_equal(result.variances.lesion_standard_deviation, np.zeros(2, dtype=np.float32))
    assert np.isfinite(result.membership.class_memberships).all()
    assert np.isfinite(result.membership.objective)
    assert np.any(result.membership.class_memberships > 0.0)
    assert result.loop_entry.gm_flag is True
    assert len(result.first_thinning.thinning.run.events) > 0
    assert result.first_thinning.thinning.run.segmentation.shape == labels.shape
    assert np.array_equal(labels, np.asanyarray(labels_nii.dataobj))

    # Reassemble the selected component calls at this seam. Individual source
    # harnesses cover those kernels; this checks that the adapter connects
    # them in the expected order and on the same constructor arrays.
    direct_relations = lesion_compute_relations(labels, 13, 3)
    np.testing.assert_array_equal(result.relations.relation_map, direct_relations.relation_map)
    np.testing.assert_array_equal(
        result.relations.relation_membership, direct_relations.relation_membership
    )
    relation_table = direct_relations.relation_membership[direct_relations.relation_map]
    variances = compute_lesion_toads_variances(
        images,
        result.initial_class_memberships,
        labels,
        class_centroids,
        result.initial_lesion_membership,
        result.initial_black_hole_membership,
        np.asarray((0.57, 0.63), dtype=np.float32),
        np.asarray((0.34, 0.4), dtype=np.float32),
        STANDARD13_TYPES,
    )
    weights = compute_lesion_toads_channel_weights(
        images,
        result.initial_class_memberships,
        class_centroids,
        result.initial_lesion_membership,
        result.initial_black_hole_membership,
        result.initial_outlier_membership,
        np.asarray((0.57, 0.63), dtype=np.float32),
        np.asarray((0.34, 0.4), dtype=np.float32),
        STANDARD13_TYPES,
        outlier_distance=0.1,
        norm_type="Identity",
    )
    np.testing.assert_array_equal(result.variances.class_standard_deviation, variances.class_standard_deviation)
    np.testing.assert_array_equal(result.variances.lesion_standard_deviation, variances.lesion_standard_deviation)
    np.testing.assert_array_equal(result.channel_weights, weights)
    ventricle_mask, ventricle = lesion_toads_ventricle_distance(labels, STANDARD13_NAMES, 2)
    gm_mask, gm = lesion_toads_gm_distance(labels, STANDARD13_NAMES, 3)
    brainstem_mask, brainstem = lesion_toads_brainstem_distance(labels, STANDARD13_NAMES, 25)
    np.testing.assert_array_equal(result.distance_maps.ventricle_mask, ventricle_mask)
    np.testing.assert_array_equal(result.distance_maps.ventricle, ventricle)
    np.testing.assert_array_equal(result.distance_maps.gray_matter_mask, gm_mask)
    np.testing.assert_array_equal(result.distance_maps.gray_matter, gm)
    np.testing.assert_array_equal(result.distance_maps.brainstem_mask, brainstem_mask)
    np.testing.assert_array_equal(result.distance_maps.brainstem, brainstem)
    expected_membership = run_lesion_toads_standard13_volume(
        images,
        class_centroids,
        variances.class_standard_deviation,
        np.asarray((0.57, 0.63), dtype=np.float32),
        variances.lesion_standard_deviation,
        weights,
        shapes,
        np.ones(13, dtype=np.float32),
        np.ones(13, dtype=np.float32),
        relation_table,
        result.initial_class_memberships,
        result.initial_lesion_membership,
        result.initial_black_hole_membership,
        result.initial_outlier_membership,
        result.initial_lesion_relation,
        ventricle,
        gm,
        brainstem,
        outlier_energy=0.1,
        lesion_power_exponent=1.2,
        gm_flag=True,
        spread=1.0,
        source_zero_variance_nan=True,
    )
    np.testing.assert_array_equal(result.membership.class_memberships, expected_membership.class_memberships)
    np.testing.assert_array_equal(result.membership.lesion_membership, expected_membership.lesion_membership)
    np.testing.assert_array_equal(result.membership.denominator, expected_membership.denominator)
    assert result.membership.objective == expected_membership.objective
    direct_thinning = run_lesion_initial_thinning_pass(
        result.loop_entry,
        labels,
        STANDARD13_TYPES,
        shapes,
        labels.copy(),
        np.zeros((*labels.shape, 13), dtype=np.bool_),
        np.zeros(labels.shape, dtype=np.float32),
        lesion_class=10,
        volume_ratio=0.5,
        limit=1000.0,
    )
    np.testing.assert_array_equal(
        result.first_thinning.thinning.run.segmentation,
        direct_thinning.thinning.run.segmentation,
    )
