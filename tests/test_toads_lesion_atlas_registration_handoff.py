"""Compose the Java-referenced 13-class pass into selected rigid registration."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_lesion_atlas_registration_handoff import (
    toads_lesion_atlas_registration_handoff,
)
from pymacruise.toads_registration_sampling import register_toads_atlas_rigid
from pymacruise.toads_shapes import _rotation_matrix, compute_toad_shape_bounds
from test_toads_lesion_atlas_type_membership import (
    LESION_TOADS_CLASS_TYPES,
    ROOT as MEMBERSHIP_FIXTURE,
    _fixture as membership_fixture,
    _parameters as membership_parameters,
    _run as java_membership_pass,
)


REFERENCE = Path(__file__).parent / "data" / "synthetic" / "toads_lesion_atlas_registration_handoff" / "java_reference.json"
JAVA = json.loads(REFERENCE.read_text())


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _from_bits(values: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32).reshape(shape)


def _args():
    data, parameters, cases = membership_fixture()
    atlas_image = nib.load(str(MEMBERSHIP_FIXTURE / "atlas_template_shapes.nii.gz"))
    subject_image = nib.load(str(MEMBERSHIP_FIXTURE / "images.nii.gz"))
    assert nib.aff2axcodes(atlas_image.affine) == nib.aff2axcodes(subject_image.affine) == ("R", "A", "S")
    atlas_shapes = np.asarray(atlas_image.dataobj, dtype=np.float32)
    selected_case = next(c for c in cases if c["field"] and c["psf"] and c["smoothing"] and c["prior"])
    membership = java_membership_pass(data, parameters, selected_case)
    options = {
        "image_resolutions_mm": tuple(float(x) for x in nib.affines.voxel_sizes(subject_image.affine)[:3]),
        "template_resolutions_mm": tuple(float(x) for x in nib.affines.voxel_sizes(atlas_image.affine)[:3]),
        "initial_transform": (0.004, -0.003, 0.002, 0.04, -0.03, 0.02),
        "offset": 0,
        "subsample": 1,
        "max_iterations": 4,
        "energy_tolerance": 0.0,
        "minimum_lambda": 1e-3,
    }
    return data, selected_case, membership, atlas_shapes, subject_image, options


def _direct_registration(memberships: np.ndarray, atlas_shapes: np.ndarray, options: dict):
    initial = np.asarray(options["initial_transform"], dtype=np.float32)
    bounds = tuple(
        compute_toad_shape_bounds(
            atlas_shapes[..., k], image_dimensions_xyz=memberships.shape[:3],
            image_resolutions_mm=options["image_resolutions_mm"],
            template_resolutions_mm=options["template_resolutions_mm"],
            transform=tuple(float(x) for x in initial),
        )
        for k in range(13)
    )
    low = np.asarray([b.lower_xyz for b in bounds], dtype=np.int32)
    high = np.asarray([b.upper_xyz for b in bounds], dtype=np.int32)
    registered = np.asarray([kind != "MASK" for kind in LESION_TOADS_CLASS_TYPES], dtype=bool)
    result = register_toads_atlas_rigid(
        memberships, atlas_shapes,
        image_resolutions_mm=options["image_resolutions_mm"],
        template_resolutions_mm=options["template_resolutions_mm"],
        initial_transform=tuple(float(x) for x in initial),
        class_bounds_low_xyz=low, class_bounds_high_xyz=high,
        registered_classes=registered,
        offset=options["offset"], subsample=options["subsample"],
        max_iterations=options["max_iterations"],
        energy_tolerance=options["energy_tolerance"],
        minimum_lambda=options["minimum_lambda"],
    )
    return result, bounds


def test_lesion_memberships_feed_registered_13_class_atlas_in_wrapper_order() -> None:
    data, case, membership, atlas_shapes, subject_image, options = _args()
    assert JAVA["wrapper_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert JAVA["atlas_source_sha256"] == "e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2"
    assert JAVA["class_types"] == list(LESION_TOADS_CLASS_TYPES)
    assert JAVA["call_order"] == [
        "computeExactMemberships", "energyRatio", "registerShapes",
        "convergenceCheck", "iterationLimitCheck",
    ]
    assert JAVA["registration_dispatch"] == "SINGLE/registerAllShapes"
    assert JAVA["registered_class_indices"] == list(range(1, 13))
    assert JAVA["registration_resets"] == ["lambda=INIT_LAMBDA", "Nturn=0", "itPlus=0"]
    assert JAVA["registration_stop"] == ["n < itMax", "diff > minEdiff", "lambda > minLambda"]
    assert JAVA["registration_finalize"] == ["computeRotation(transform)", "computeTransformedShapeBoundingBox()"]

    # The 16-case method reference is Java-produced by the P18czl harness. Feed
    # those exact membership bits to the independent Python composition so this
    # check does not merely register the port's own pass output.
    assert _bits(membership.memberships) == case["membership_bits"]
    java_memberships = _from_bits(case["membership_bits"], membership.memberships.shape)
    expected, initial_bounds = _direct_registration(java_memberships, atlas_shapes, options)
    before_memberships = membership.memberships.copy()
    before_lesion_map = membership.lesion_membership.copy()
    result = toads_lesion_atlas_registration_handoff(
        membership, atlas_shapes, **options,
    )
    assert result.registration is not None
    assert result.registered_class_indices == tuple(range(1, 13))
    np.testing.assert_array_equal(result.membership_pass.memberships, before_memberships)
    np.testing.assert_array_equal(result.membership_pass.lesion_membership, before_lesion_map)
    np.testing.assert_array_equal(result.registration.transform, expected.transform)
    assert _bits(result.registration.energy) == _bits(expected.energy)
    assert result.registration.stop_reason == expected.stop_reason == "iteration_limit"
    assert len(result.registration.iterations) == len(expected.iterations) == options["max_iterations"]
    for actual, reference in zip(result.registration.iterations, expected.iterations):
        np.testing.assert_array_equal(actual.transform, reference.transform)
        assert _bits(actual.energy) == _bits(reference.energy)
        assert _bits(actual.lambda_value) == _bits(reference.lambda_value)
        assert _bits(actual.relative_change) == _bits(reference.relative_change)
        assert (actual.accepted, actual.rejected_trials) == (reference.accepted, reference.rejected_trials)
    assert len(result.shape_bounds) == 13
    assert result.shape_bounds == tuple(
        compute_toad_shape_bounds(
            atlas_shapes[..., k], image_dimensions_xyz=membership.memberships.shape[:3],
            image_resolutions_mm=options["image_resolutions_mm"],
            template_resolutions_mm=options["template_resolutions_mm"],
            transform=tuple(float(x) for x in result.transform),
        )
        for k in range(13)
    )
    assert any(b.support_count > 0 for b in initial_bounds[10:11])
    assert result.rotation.shape == (3, 3)
    np.testing.assert_allclose(result.rotation, _rotation_matrix(result.transform[:3]), rtol=0, atol=0)
    # Lesion slot 10 must participate: removing only its freshly updated
    # membership changes the exact registration state on this folded fixture.
    without_lesion = java_memberships.copy()
    without_lesion[..., 10] = 0
    no_lesion_result, _ = _direct_registration(without_lesion, atlas_shapes, options)
    assert not np.array_equal(result.transform, no_lesion_result.transform)
    out_image = nib.Nifti1Image(result.membership_pass.memberships, subject_image.affine)
    assert nib.aff2axcodes(out_image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(out_image.affine, subject_image.affine)


def test_lesion_atlas_alignment_disabled_keeps_transform_and_class_contract() -> None:
    _, _, membership, atlas_shapes, _, options = _args()
    result = toads_lesion_atlas_registration_handoff(
        membership, atlas_shapes, align_shapes=False, **options,
    )
    assert result.registration is None
    assert result.registered_class_indices == tuple(range(1, 13))
    np.testing.assert_array_equal(result.transform, np.asarray(options["initial_transform"], dtype=np.float32))
    np.testing.assert_allclose(result.rotation @ result.rotation.T, np.eye(3), atol=2e-6)
    assert len(result.shape_bounds) == 13


@pytest.mark.parametrize("bad", [None, object()])
def test_registration_handoff_rejects_non_membership_results(bad) -> None:
    _, _, _, atlas_shapes, _, options = _args()
    with pytest.raises(TypeError, match="ToadsLesionAtlasTypeMembershipPass"):
        toads_lesion_atlas_registration_handoff(bad, atlas_shapes, **options)


def test_registration_handoff_rejects_template_class_count_mismatch() -> None:
    _, _, membership, atlas_shapes, _, options = _args()
    with pytest.raises(ValueError, match="thirteen class maps"):
        toads_lesion_atlas_registration_handoff(membership, atlas_shapes[..., :-1], **options)
