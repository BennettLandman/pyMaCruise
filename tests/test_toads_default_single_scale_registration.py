"""Source-pinned later atlas-aware membership and scale-one registration slice."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from pymacruise.toads_default_single_scale_registration import (
    initialize_and_register_toads_default_single_scale,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_single_scale_registration"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsDefaultSingleScaleRegistrationReferenceHarness.java"
HARNESS_SHA256 = "32f99e62674d2712829cd6ac1d57d5f4bc196368c8aef625cce73e5bb15889bf"
SOURCE_HASHES = {
    "wrapper_source_sha256": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
    "preprocess_source_sha256": "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7",
    "atlas_source_sha256": "e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2",
}
INITIAL = (.031, -.027, .022, .31, -.27, .22)
RESOLUTIONS = (.9, 1.2, 1.45)


def _bits(value: np.ndarray | np.float32) -> list[int]:
    return np.asarray(value, dtype=np.float32).ravel().view(np.int32).tolist()


def _decode_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32)


def _ordered_bits(values: np.ndarray | np.float32) -> np.ndarray:
    signed = np.asarray(values, dtype=np.float32).ravel().view(np.int32).astype(np.int64)
    return np.where(signed < 0, -signed, signed + (1 << 31))


def _max_ulp_difference(actual: np.ndarray | np.float32, expected_bits: list[int]) -> int:
    return int(np.max(np.abs(_ordered_bits(actual) - _ordered_bits(_decode_bits(expected_bits)))))


def _inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    volumes = tuple(nib.load(str(ROOT / f"{name}.nii.gz")) for name in (
        "images", "transformed_shapes", "atlas_shapes",
    ))
    for volume in volumes:
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
        np.testing.assert_allclose(nib.affines.voxel_sizes(volume.affine), RESOLUTIONS, rtol=0, atol=1e-6)
    assert volumes[0].shape == (9, 8, 7, 2)
    assert volumes[1].shape == (9, 8, 7, 11)
    assert volumes[2].shape == (11, 10, 9, 11)
    return tuple(np.asarray(volume.dataobj, dtype=np.float32) for volume in volumes)


def _centroids() -> np.ndarray:
    return np.asarray(
        [[np.float32(.19) + np.float32(.021) * k + np.float32(.013) * channel
          for k in range(11)] for channel in range(2)], dtype=np.float32,
    )


def _run(images: np.ndarray, transformed_shapes: np.ndarray, atlas_shapes: np.ndarray, centroids: np.ndarray,
         *, align_shapes: bool = True):
    return initialize_and_register_toads_default_single_scale(
        images, centroids, transformed_shapes, atlas_shapes, approx=JAVA["approx"],
        initial_transform=INITIAL, image_resolutions_mm=RESOLUTIONS,
        template_resolutions_mm=RESOLUTIONS, max_iterations=1, align_shapes=align_shapes,
    )


def test_atlas_aware_single_scale_registration_matches_pinned_java_slice() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    for key, expected_sha in SOURCE_HASHES.items():
        assert JAVA[key] == expected_sha
    assert JAVA["class_types"] == list(DEFAULT_ATLAS_CLASS_TYPES)
    assert JAVA["call_order"] == [
        "initialApproxMemberships(atlas)",
        "updateShapeRegistration(rigid,iter=1,level=1)",
        "registerShapes()->registerAllShapes(single-scale)",
    ]

    images, transformed_shapes, atlas_shapes = _inputs()
    inputs_before = tuple(value.copy() for value in (images, transformed_shapes, atlas_shapes))
    result = _run(images, transformed_shapes, atlas_shapes, _centroids())
    assert result.call_order == (
        "preprocess.initialMemberships(memberships,best,atlas,approx,centroids,classes)",
        'atlas.updateShapeRegistration(memberships,best,segmentation,modality,"rigid",alignIter,1,...)',
        "atlas.registerShapes() -> registerAllShapes()",
    )
    assert _bits(result.initialization.memberships) == JAVA["initial_memberships_bits"]
    np.testing.assert_array_equal(
        result.initialization.best_class_indices.ravel().astype(int), JAVA["best_indices"],
    )
    assert np.count_nonzero(result.initialization.best_class_indices == 0) > 0
    assert result.registration is not None
    assert len(result.registration.iterations) == len(JAVA["iterations"]) == 1
    iteration = result.registration.iterations[0]
    reference = JAVA["iterations"][0]
    gradient, hessian = result.iteration_coefficients[0]
    # Measured on this scale-one fixture: the only coefficient drift is at
    # most 8 ULP for gradients and 2 ULP for Hessians.
    assert _max_ulp_difference(gradient, reference["gradient"]) <= 8
    assert _max_ulp_difference(hessian, reference["hessian"]) <= 2
    trials = result.iteration_trials[0]
    assert len(trials) == len(reference["trials"]) == 1
    trial_transform, trial_energy = trials[0]
    java_trial = reference["trials"][0]
    assert _max_ulp_difference(trial_transform, java_trial["transform"]) <= 16
    assert _max_ulp_difference(trial_energy, [java_trial["energy"]]) <= 2
    assert iteration.accepted is reference["accepted"] is True
    assert iteration.rejected_trials == reference["rejected_trials"] == 0
    assert _bits(iteration.lambda_value) == [reference["lambda"]]
    assert _max_ulp_difference(iteration.energy, [reference["energy"]]) <= 2
    assert abs(float(iteration.relative_change) - float(_decode_bits([reference["relative_change"]])[0])) <= 3.3e-7
    assert _max_ulp_difference(result.transform, JAVA["final_transform_bits"]) <= 16
    np.testing.assert_allclose(
        result.rotation, _decode_bits(JAVA["rotation_bits"]).reshape(3, 3),
        rtol=0, atol=6.3e-8,
    )
    py_bounds = [
        [*bound.lower_xyz, *bound.upper_xyz, bound.support_count]
        for bound in result.shape_bounds
    ]
    assert py_bounds == JAVA["bounds"]
    assert _bits(iteration.transform) == _bits(trial_transform)
    for actual, before in zip((images, transformed_shapes, atlas_shapes), inputs_before, strict=True):
        np.testing.assert_array_equal(actual, before)


def test_update_without_optional_shape_alignment_and_candidates_affect_registration() -> None:
    images, transformed_shapes, atlas_shapes = _inputs()
    centroids = _centroids()
    originals = tuple(value.copy() for value in (images, transformed_shapes, atlas_shapes))
    no_align = _run(images, transformed_shapes, atlas_shapes, centroids, align_shapes=False)
    assert no_align.registration is None
    np.testing.assert_array_equal(no_align.transform, np.asarray(INITIAL, dtype=np.float32))
    assert len(no_align.call_order) == 2
    assert not no_align.iteration_coefficients and not no_align.iteration_trials

    changed_centroids = centroids + np.float32(.05)
    changed = _run(images, transformed_shapes, atlas_shapes, changed_centroids)
    assert np.count_nonzero(
        no_align.initialization.best_class_indices != changed.initialization.best_class_indices
    ) > 0
    aligned = _run(images, transformed_shapes, atlas_shapes, centroids)
    assert float(np.max(np.abs(aligned.transform - changed.transform))) > .1
    for actual, before in zip((images, transformed_shapes, atlas_shapes), originals, strict=True):
        np.testing.assert_array_equal(actual, before)
