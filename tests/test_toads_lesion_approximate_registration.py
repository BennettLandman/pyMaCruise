"""Selected-Java comparison for top-k initialization and the approximate pyramid."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_approximate_registration import (
    initialize_toads_approximate_memberships,
    register_toads_approximate_lesion_atlas_pyramid,
    subsample_toads_approximate_memberships,
)
from pymacruise.toads_lesion_atlas_registration_handoff import (
    toads_lesion_atlas_registration_handoff,
)
from test_toads_lesion_atlas_type_membership import (
    LESION_TOADS_CLASS_TYPES,
    ROOT as MEMBERSHIP_FIXTURE,
    _fixture,
    _parameters,
    _run,
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


def test_selected_approximate_memberships_and_lesion_atlas_pyramid_match_java() -> None:
    data, params, cases, arrays, subject, template = _inputs()
    assert JAVA["wrapper_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert JAVA["preprocess_source_sha256"] == "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7"
    assert JAVA["atlas_source_sha256"] == "e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2"
    assert JAVA["class_types"] == list(LESION_TOADS_CLASS_TYPES)
    assert JAVA["approx"] == 3 and JAVA["levels"] == 3 and JAVA["subsample"] == 2
    assert JAVA["call_order"] == [
        "initialMemberships", "initShapeRegistration", "setTransform(preprocess)",
        "updateShapeRegistration", "registerShapesPyramid",
    ]
    assert JAVA["pyramid_dispatch"] == "registerShapesApproxPyramid"
    assert JAVA["scale_order"] == [4, 2]
    assert JAVA["registered_class_indices"] == list(range(1, 13))
    assert DESCRIPTOR.is_file()

    before = {key: value.copy() for key, value in arrays.items()}
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
    omitted_with_subject_and_template_support = False
    for point in np.ndindex(initialization.memberships.shape[:3]):
        for klass in range(13):
            if klass not in initialization.best_class_indices[point] and arrays["transformed_shapes"][point + (klass,)] > 0 and np.any(arrays["atlas_shapes"][..., klass] > 0):
                omitted_with_subject_and_template_support = True
                break
        if omitted_with_subject_and_template_support:
            break
    assert omitted_with_subject_and_template_support

    result = register_toads_approximate_lesion_atlas_pyramid(
        initialization,
        arrays["atlas_shapes"],
        image_resolutions_mm=(.9, 1.2, 1.45),
        template_resolutions_mm=(.9, 1.2, 1.45),
        initial_transform=(.004, -.003, .002, .04, -.03, .02),
        levels=3,
        max_iterations=3,
    )
    assert len(result.levels) == 2
    assert len(JAVA["stages"]) == 2
    trial_energy_ulp_differences: list[int] = []
    for actual, reference in zip(result.levels, JAVA["stages"]):
        assert actual.scale == reference["scale"]
        assert list(actual.memberships.shape[:3]) == reference["dims"]
        assert list(actual.result.iterations[0].transform.shape) == [6]
        assert _bits(actual.memberships) == reference["memberships_bits"]
        np.testing.assert_array_equal(actual.best_class_indices.ravel().astype(int), reference["best_indices"])
        assert _bits(actual.initial_transform) == reference["initial_transform_bits"]
        assert len(actual.iteration_coefficients) == len(reference["iterations"])
        assert len(actual.iteration_trials) == len(reference["iterations"])
        for coeff, trials, expected_iteration, actual_iteration in zip(
            actual.iteration_coefficients, actual.iteration_trials,
            reference["iterations"], actual.result.iterations,
        ):
            assert _bits(coeff[0]) == expected_iteration["gradient"]
            assert _bits(coeff[1]) == expected_iteration["hessian"]
            assert len(trials) == len(expected_iteration["trials"])
            for (trial_transform, trial_energy), expected_trial in zip(trials, expected_iteration["trials"]):
                assert _bits(trial_transform) == expected_trial["transform"]
                # Image interpolation and accumulation can differ by one
                # float32 rounding step across Java and NumPy. Keep the
                # allowance explicit and tiny; the exact source line-search
                # accept/reject decision and resulting state are asserted
                # below, so a branch-changing discrepancy still fails.
                actual_energy_bits = _bits(trial_energy)[0]
                trial_energy_ulp_differences.append(
                    abs(actual_energy_bits - expected_trial["energy"])
                )
            assert _bits(actual_iteration.transform) == expected_iteration["transform"]
            assert _bits(actual_iteration.energy) == [expected_iteration["energy"]]
            assert _bits(actual_iteration.lambda_value) == [expected_iteration["lambda"]]
            assert _bits(actual_iteration.relative_change) == [expected_iteration["relative_change"]]
            assert actual_iteration.accepted == expected_iteration["accepted"]
            assert actual_iteration.rejected_trials == expected_iteration["rejected_trials"]
        assert _bits(actual.result.transform) == reference["transform_bits"]
        assert _bits(actual.result.energy) == [reference["energy_bits"]]
        assert _bits(actual.result.iterations[-1].lambda_value) == [reference["lambda_bits"]]
        assert _bits(actual.result.iterations[-1].relative_change) == [reference["relative_change_bits"]]
    assert _bits(result.transform) == JAVA["final_transform_bits"]
    assert trial_energy_ulp_differences
    assert max(trial_energy_ulp_differences) <= 1
    assert len(result.shape_bounds) == 13
    assert result.shape_bounds[10].support_count > 0
    assert result.rotation.shape == (3, 3)
    nifti = nib.Nifti1Image(arrays["images"], subject.affine)
    assert nib.aff2axcodes(nifti.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(nifti.affine, subject.affine)
    for key, value in before.items():
        np.testing.assert_array_equal(arrays[key], value)

    # The selected exact membership pass provides a separate, later starting
    # route (P18czm); top-k initialization must produce a materially distinct
    # registration trajectory on this fixture.
    selected_case = next(c for c in cases if c["field"] and c["psf"] and c["smoothing"] and c["prior"])
    exact_membership = _run(data, params, selected_case)
    exact_registration = toads_lesion_atlas_registration_handoff(
        exact_membership, arrays["atlas_shapes"],
        image_resolutions_mm=(.9, 1.2, 1.45),
        template_resolutions_mm=(.9, 1.2, 1.45),
        initial_transform=(.004, -.003, .002, .04, -.03, .02),
        max_iterations=3,
        energy_tolerance=0.0,
    )
    assert not np.array_equal(result.transform, exact_registration.transform)


def test_ties_force_shape_winner_then_lowest_remaining_class_indices() -> None:
    image = np.zeros((3, 3, 3, 1), dtype=np.float32)
    centroids = np.zeros((1, 13), dtype=np.float32)
    priors = np.zeros((3, 3, 3, 13), dtype=np.float32)
    priors[..., 5] = 0.9
    result = initialize_toads_approximate_memberships(
        image, centroids, priors, LESION_TOADS_CLASS_TYPES, approx=3,
    )
    np.testing.assert_array_equal(result.best_class_indices[1, 1, 1], [5, 0, 1])


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
