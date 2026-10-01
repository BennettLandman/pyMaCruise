from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_initial_registration_control import register_initial_memberships_and_control
from pymacruise.toads_initial_thinning_composition import toads_initial_fielded_membership_to_thinning
from pymacruise.toads_initial_thinning_statistics_membership import continue_initial_thinning_through_memberships


ROOT = Path(__file__).parent / "data/synthetic/toads_initial_thinning_startup"
JAVA_P18CYO_REFERENCE = ROOT / "continuation_java_reference.json"


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _p18cyo():
    spec = json.loads((ROOT / "spec.json").read_text())
    p = spec["parameters"]
    image = _load("images_ras.nii.gz").astype(np.float32)
    segmentation = _load("segmentation_ras.nii.gz").astype(np.uint8)
    memberships = _load("initial_memberships_ras.nii.gz").astype(np.float32)
    fields = _load("fields_ras.nii.gz").astype(np.float32)
    bias = _load("psf_bias_ras.nii.gz").astype(np.float32)
    transformed = _load("transformed_shapes_ras.nii.gz").astype(np.float32)
    centers = np.asarray(p["centroids"], dtype=np.float32)
    exponents = np.asarray(p["exponent_by_channel_class"], dtype=np.float32)
    startup = toads_initial_fielded_membership_to_thinning(
        image_values=image, segmentation=segmentation, initial_memberships=memberships,
        centroids=centers, exponent_by_channel_class=exponents,
        modality_weights=p["modality_weights"], class_variance_weights=p["class_variance_weights"],
        image_variances=p["image_variances"], transformed_shapes=transformed,
        smoothing_priors=p["smoothing_priors"], centroid_similarity=p["centroid_similarity"],
        atlas_shape_priors=p["atlas_shape_priors"], smoothing=spec["membership_smoothing"],
        prior_factor=spec["prior_factor"], prior_coefficient=spec["prior_coefficient"],
        channel_fields=[fields[..., 0], fields[..., 1]], limit=spec["relation_limit"],
        volume_ratio=spec["thinning_volume_ratio"], thinning_limit=spec["thinning_limit"], psf_bias=bias,
    )
    category_exponents = np.asarray([exponents[:, 1], exponents[:, 3], exponents[:, 8]], dtype=np.float32).T
    continued = continue_initial_thinning_through_memberships(
        startup, image_values=image, initial_centroids=centers,
        category_exponents=category_exponents, transformed_shapes=transformed,
        image_variance_weights=p["image_variances"], smoothing_priors=p["smoothing_priors"],
        atlas_shape_priors=p["atlas_shape_priors"], exponent_by_channel_class=exponents,
        modality_weights=p["modality_weights"], prior_factor=spec["prior_factor"],
        prior_coefficient=spec["prior_coefficient"], smoothing=spec["membership_smoothing"],
        prior_scale=0.05, intensity_field=fields, psf_bias=bias,
    )
    return continued, transformed


def test_p18cyo_membership_and_objective_feed_registration_then_control() -> None:
    continuation, atlas = _p18cyo()
    source_reference = json.loads(JAVA_P18CYO_REFERENCE.read_text())
    membership_digest = hashlib.sha256(
        np.asarray(continuation.membership_pass.memberships, dtype=">f4").tobytes(order="C")
    ).hexdigest()
    assert membership_digest == source_reference["updated_membership_sha256"]
    assert int(np.asarray(continuation.membership_pass.objective, dtype=np.float32).view(np.int32)) == source_reference["objective_bits"]
    membership_before = continuation.membership_pass.memberships.copy()
    objective = np.float32(continuation.membership_pass.objective)
    result = register_initial_memberships_and_control(
        continuation, atlas, image_resolutions_mm=(1, 1, 1), template_resolutions_mm=(1, 1, 1),
        initial_transform=(0, 0, 0, 0, 0, 0), previous_energy=np.float32(objective * 1.01),
        iteration=3, maximum_iterations=5, convergence_threshold=1e-4, do_growing=False,
    )
    assert result.handoff.membership_pass is continuation.membership_pass
    np.testing.assert_array_equal(result.handoff.membership_pass.memberships, membership_before)
    assert result.handoff.registration is not None
    assert len(result.handoff.registration.iterations) in range(1, 6)
    assert result.handoff.transform.shape == (6,)
    assert result.handoff.rotation.shape == (3, 3)
    assert len(result.handoff.shape_bounds) == 11
    assert result.control.energy_ratio == np.float32(
        np.abs(np.float32(np.float32(2) * np.float32(np.float32(np.float32(objective * 1.01) - objective)) /
                          np.float32(np.float32(objective * 1.01) + objective)))
    )
    assert result.stage_order[-4:] == (
        "atlas.registerShapes() -> registerAllShapes()", "energyRatio", "convergenceCheck", "iterationLimitCheck",
    )
