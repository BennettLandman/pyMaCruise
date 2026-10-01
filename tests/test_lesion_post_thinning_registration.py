from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_energy_ratio import lesion_toads_energy_ratio
from pymacruise.lesion_post_thinning_registration import register_after_post_thinning_membership
from test_lesion_post_thinning_gm_membership import (
    _bits,
    _compose,
    SPEC as UPSTREAM_SPEC,
)


ROOT = Path(__file__).parent / "data/synthetic/lesion_post_thinning_registration"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))
SPEC = json.loads((ROOT / "fixture.json").read_text(encoding="utf-8"))
HARNESS = Path(__file__).parent / "java_reference/LesionPostThinningRegistrationReferenceHarness.java"
MEMBERSHIP_REFERENCE = json.loads(
    (Path(__file__).parent / "data/synthetic/lesion_post_thinning_gm_membership_update/java_reference.json").read_text()
)
RATIO_REFERENCE = json.loads(
    (Path(__file__).parent / "data/synthetic/lesion_energy_ratio/java_reference.json").read_text()
)


def test_post_thinning_registration_matches_java_trace_even_when_ratio_is_nan() -> None:
    # This reproduces P39ak→P39am in source-faithful mode. P39an's energy ratio
    # is NaN, but the selected wrapper calls registerShapes before convergence.
    from pymacruise.lesion_post_thinning_membership import run_post_thinning_standard13_membership_update

    frames, state, initial, thinning = _compose()
    membership_update = run_post_thinning_standard13_membership_update(
        thinning,
        state["images"],
        UPSTREAM_SPEC["class_centroids"],
        state["lesion"], state["black"], state["outlier"],
        UPSTREAM_SPEC["lesion_centroid"], (0.0, 0.0),
        state["shape_values"], UPSTREAM_SPEC["shape_priors"], UPSTREAM_SPEC["smoothing_priors"],
        state["direct_relations"], frames[..., 63], state["ventricle_distances"],
        state["gray_matter_distances"], state["brainstem_distances"], UPSTREAM_SPEC["object_types"],
        centroid_options={"mode": "none", "modalities": ("T1_MPRAGE", "FLAIR"),
                          "lesion_power_exponents": (0.4, 1.2)},
        channel_weight_options={"outlier_distance": UPSTREAM_SPEC["outlier_energy"]},
        outlier_energy=UPSTREAM_SPEC["outlier_energy"],
        lesion_power_exponent=UPSTREAM_SPEC["lesion_power_exponent"], flair_channel=1,
    )
    assert _bits(membership_update.membership.class_memberships) == MEMBERSHIP_REFERENCE["class_membership_bits"]
    assert _bits(membership_update.membership.objective) == MEMBERSHIP_REFERENCE["objective_bits"]
    assert MEMBERSHIP_REFERENCE["objective_bits"] == 0x7FC00000

    ratio = lesion_toads_energy_ratio(initial.objective, membership_update.membership.objective)
    ratio_case = next(case for case in RATIO_REFERENCE["cases"] if case["case"] == "p39ak_to_p39am_nan")
    assert f"{int(np.asarray(ratio.energy_ratio, dtype=np.float32).view(np.uint32)):08x}" == ratio_case["energy_ratio_bits"]
    assert ratio_case["energy_ratio_bits"] == "ffc00000"

    image = nib.load(str(ROOT / "registration_inputs_ras.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert list(image.shape) == SPEC["dimensions"]
    frames26 = np.asarray(image.dataobj, dtype=np.float32)
    memberships = frames26[..., :13]
    shapes = tuple(frames26[..., 13 + index] for index in range(13))
    np.testing.assert_array_equal(memberships, membership_update.membership.class_memberships)
    assert np.isfinite(memberships).all()
    assert hashlib.sha256((ROOT / "registration_inputs_ras.nii.gz").read_bytes()).hexdigest() == SPEC["fixture_sha256"]
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == REFERENCE["harness_sha256"]

    result = register_after_post_thinning_membership(
        membership_update,
        ratio,
        shapes,
        SPEC["initial_transform"],
        image_resolutions_mm=tuple(SPEC["image_resolutions_mm"]),
        template_resolutions_mm=tuple(SPEC["template_resolutions_mm"]),
        align_shapes=True,
        init_align_iter=SPEC["init_align_iter"],
    )
    assert result.stage_order[-2:] == ("energyRatio", "registerShapes")
    assert np.isnan(result.energy_ratio.energy_ratio)
    optimization = result.registration.optimization
    assert optimization is not None, "NaN ratio must not suppress registerShapes()"
    assert optimization.stop_reason == REFERENCE["stop_reason"]
    assert len(optimization.iterations) == len(REFERENCE["iterations"])
    np.testing.assert_array_equal(result.registration.transform, np.asarray(REFERENCE["transform"], dtype=np.float32))
    _assert_f32(optimization.energy, REFERENCE["energy"], "final energy")
    _assert_f32(optimization.lambda_value, REFERENCE["lambda"], "final lambda")

    for step, expected in zip(optimization.iterations, REFERENCE["iterations"]):
        assert len(step.trial_energies) == len(expected["trials"])
        assert step.result.accepted == expected["accepted"]
        assert step.result.rejected_trials == expected["rejected_trials"]
        np.testing.assert_array_equal(step.result.transform, np.asarray(expected["transform"], dtype=np.float32))
        _assert_f32(step.result.energy, expected["energy"], "iteration energy")
        _assert_f32(step.result.relative_change, expected["relative_change"], "relative change")
        _assert_f32(step.result.lambda_value, expected["lambda"], "iteration lambda")
        for trial_index, expected_trial in enumerate(expected["trials"]):
            np.testing.assert_array_equal(
                step.trial_transforms[trial_index], np.asarray(expected_trial["transform"], dtype=np.float32)
            )
            _assert_f32(step.trial_energies[trial_index], expected_trial["energy"], "trial energy")

    np.testing.assert_array_equal(result.registration.rotation, np.asarray(REFERENCE["rotation"], dtype=np.float32))
    for bounds, expected in zip(result.registration.shape_bounds, REFERENCE["bounds"]):
        assert bounds.lower_xyz == tuple(expected["lower"])
        assert bounds.upper_xyz == tuple(expected["upper"])
        assert bounds.support_count == expected["support_count"]


def _assert_f32(actual: float, expected: float, label: str) -> None:
    actual_bits = int(np.asarray(actual, dtype=np.float32).view(np.uint32))
    expected_bits = int(np.asarray(expected, dtype=np.float32).view(np.uint32))
    assert actual_bits == expected_bits, label
