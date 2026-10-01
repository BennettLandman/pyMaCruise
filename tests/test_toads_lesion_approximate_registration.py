"""Selected-Java validation of the atlas-aware lesion top-k overload."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_approximate_registration import (
    initialize_toads_approximate_memberships,
    subsample_toads_approximate_memberships,
)
from test_toads_lesion_atlas_type_membership import (
    LESION_TOADS_CLASS_TYPES,
    ROOT as MEMBERSHIP_FIXTURE,
    _fixture,
    _parameters,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_lesion_approximate_registration"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
DESCRIPTOR = Path(__file__).parents[1] / "data" / "atlases" / "lesionToads-atlas-2012" / "cruise-atlas-12obj-lesiontoads2012.txt"


def _bits(values: np.ndarray | np.float32) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _inputs():
    data, params, cases = _fixture()
    subject = nib.load(str(MEMBERSHIP_FIXTURE / "images.nii.gz"))
    transformed = nib.load(str(MEMBERSHIP_FIXTURE / "transformed_shapes.nii.gz"))
    template = nib.load(str(MEMBERSHIP_FIXTURE / "atlas_template_shapes.nii.gz"))
    assert nib.aff2axcodes(subject.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(transformed.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(template.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(subject.affine), (.9, 1.2, 1.45), rtol=0, atol=1e-6)
    np.testing.assert_allclose(nib.affines.voxel_sizes(template.affine), (.9, 1.2, 1.45), rtol=0, atol=1e-6)
    np.testing.assert_array_equal(subject.affine, transformed.affine)
    assert subject.shape[:3] == (9, 8, 7)
    assert template.shape == (11, 10, 9, 13)
    data_arrays = {
        "images": np.asarray(subject.dataobj, dtype=np.float32),
        "transformed_shapes": np.asarray(transformed.dataobj, dtype=np.float32),
        "atlas_shapes": np.asarray(template.dataobj, dtype=np.float32),
    }
    return data, params, cases, data_arrays, subject, template


def _initial(data_arrays, params):
    return initialize_toads_approximate_memberships(
        data_arrays["images"], params["centroids"], data_arrays["transformed_shapes"],
        LESION_TOADS_CLASS_TYPES, approx=3,
    )


def test_atlas_aware_lesion_approximate_memberships_match_java() -> None:
    _, params, _, arrays, subject, template = _inputs()
    assert JAVA["wrapper_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert JAVA["preprocess_source_sha256"] == "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7"
    assert JAVA["atlas_source_sha256"] == "e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2"
    assert JAVA["lesion_wrapper_source_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    assert JAVA["class_types"] == list(LESION_TOADS_CLASS_TYPES)
    assert JAVA["approx"] == 3
    assert JAVA["validation_scope"] == "method-only atlas-aware 13-class approximate-membership overload"
    assert JAVA["selected_overload"] == (
        "preprocess.initialMemberships(segmentation.getMemberships(),segmentation.getBest(),"
        "atlas,approx,centroids,classes)"
    )
    assert JAVA["selected_wrapper"] == "AlgorithmMultipleToads"
    assert JAVA["source_call_order"] == [
        "new MultipleToadSegmentation(...)",
        "segmentation.initCentroids(centroids)",
        "segmentation.setLesionCentroid(...)",
        "preprocess.initialMemberships(segmentation.getMemberships(),segmentation.getBest(),atlas,approx,centroids,classes)",
        "atlas.updateShapeRegistration(...)",
        "atlas.registerShapes() when alignShapes",
    ]
    assert JAVA["reachability_classification"] == {
        "selected": "atlas-aware overload is reached from AlgorithmMultipleToads after segmentation initialization",
        "initial_pyramid": "earlier AlgorithmMultipleToads call is the no-atlas overload",
        "dedicated_lesion_wrapper": "AlgorithmLesionToads uses the no-atlas overload for its initial registration pyramid",
        "scope": "overload method only; no registration result claimed",
    }
    evidence = JAVA["p18czn_evidence_sha256"]
    assert evidence["harness"] == hashlib.sha256(
        (Path(__file__).parent / "java_reference/ToadsLesionApproximateRegistrationReferenceHarness.java").read_bytes()
    ).hexdigest()
    assert evidence["numerics"] == hashlib.sha256(
        (Path(__file__).parent / "java_reference/upstream/edu/jhmi/rad/medic/utilities/Numerics.java").read_bytes()
    ).hexdigest()
    for key, path in {
        "images": MEMBERSHIP_FIXTURE / "images.nii.gz",
        "transformed_shapes": MEMBERSHIP_FIXTURE / "transformed_shapes.nii.gz",
        "template_shapes": MEMBERSHIP_FIXTURE / "atlas_template_shapes.nii.gz",
        "descriptor": DESCRIPTOR,
    }.items():
        assert evidence[key] == hashlib.sha256(path.read_bytes()).hexdigest(), key
    assert DESCRIPTOR.is_file()

    arrays_before = {key: value.copy() for key, value in arrays.items()}
    params_before = {key: value.copy() for key, value in params.items() if isinstance(value, np.ndarray)}
    initialization = _initial(arrays, params)
    assert initialization.class_types == LESION_TOADS_CLASS_TYPES
    assert initialization.memberships.shape == (9, 8, 7, 3)
    assert initialization.best_class_indices.shape == (9, 8, 7, 3)
    assert _bits(initialization.memberships) == JAVA["initial_memberships_bits"]
    np.testing.assert_array_equal(
        initialization.best_class_indices.ravel().astype(int), JAVA["initial_best_indices"]
    )
    assert np.all(initialization.best_class_indices < 13)
    assert np.any(initialization.best_class_indices == 10)
    assert np.all(np.diff(np.sort(initialization.best_class_indices, axis=-1), axis=-1) > 0)
    # All provided grid locations, including every edge, are compared above
    # against the harness output; the outer input boundary is nontrivial.
    assert np.any(initialization.memberships[0, :, :, :] != 0.0)
    for key, value in arrays_before.items():
        np.testing.assert_array_equal(arrays[key], value)
    for key, value in params_before.items():
        np.testing.assert_array_equal(params[key], value)
    nifti = nib.Nifti1Image(arrays["images"], subject.affine)
    assert nib.aff2axcodes(nifti.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(nifti.affine, subject.affine)


def test_ties_force_shape_winner_then_lowest_remaining_class_indices() -> None:
    image = np.zeros((1, 1, 1, 1), dtype=np.float32)
    centroids = np.zeros((1, 13), dtype=np.float32)
    centroids[0, 5] = 1.0  # equalizes every class distance with shape[5] = 1
    priors = np.zeros((1, 1, 1, 13), dtype=np.float32)
    priors[..., 5] = 1.0
    result = initialize_toads_approximate_memberships(
        image, centroids, priors, LESION_TOADS_CLASS_TYPES, approx=3,
    )
    assert result.best_class_indices[0, 0, 0].tolist() == JAVA["shape_seeded_membership_tie"]["indices"]
    assert _bits(result.memberships[0, 0, 0]) == JAVA["shape_seeded_membership_tie"]["membership_bits"]


def test_zero_cost_reciprocal_fallback_matches_java_topk() -> None:
    result = initialize_toads_approximate_memberships(
        np.zeros((1, 1, 1, 1), dtype=np.float32),
        np.zeros((1, 13), dtype=np.float32),
        np.zeros((1, 1, 1, 13), dtype=np.float32),
        LESION_TOADS_CLASS_TYPES,
        approx=3,
    )
    assert result.best_class_indices[0, 0, 0].tolist() == JAVA["zero_cost_tie"]["indices"]
    assert _bits(result.memberships[0, 0, 0]) == JAVA["zero_cost_tie"]["membership_bits"]


def test_approximate_membership_pyramid_pooling_discards_partial_high_edges() -> None:
    mem = np.zeros((5, 4, 3, 3), dtype=np.float32)
    best = np.zeros_like(mem, dtype=np.uint8)
    mem[..., 0] = .1
    mem[..., 1] = .2
    mem[..., 2] = .3
    best[..., 0] = 12
    best[..., 1] = 2
    best[..., 2] = 10
    result = subsample_toads_approximate_memberships(mem, best, classes=13, factor=2)
    assert result.memberships.shape == (2, 2, 1, 3)
    np.testing.assert_array_equal(result.best_class_indices[0, 0, 0], [10, 2, 12])
    np.testing.assert_allclose(result.memberships[0, 0, 0], [.3, .2, .1], rtol=0, atol=1e-7)


@pytest.mark.parametrize("approx", [0, 1, 13, 14])
def test_approx_initialization_rejects_non_topk_values(approx: int) -> None:
    _, params, _, arrays, _, _ = _inputs()
    with pytest.raises(ValueError, match="approx"):
        initialize_toads_approximate_memberships(
            arrays["images"], params["centroids"], arrays["transformed_shapes"],
            LESION_TOADS_CLASS_TYPES, approx=approx,
        )
