from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_energy_ratio import lesion_toads_energy_ratio
from pymacruise.lesion_post_registration import lesion_toads_post_registration_branch
from pymacruise.lesion_post_thinning_membership import run_post_thinning_standard13_membership_update
from pymacruise.lesion_post_thinning_registration import register_after_post_thinning_membership
from test_lesion_post_thinning_gm_membership import _compose, SPEC


ROOT = Path(__file__).parent / "data/synthetic/lesion_post_registration_nan"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))
HARNESS = Path(__file__).parent / "java_reference/LesionNaNPostRegistrationBranchReferenceHarness.java"
REGISTRATION_ROOT = Path(__file__).parent / "data/synthetic/lesion_post_thinning_registration"
REGISTRATION_SPEC = json.loads((REGISTRATION_ROOT / "fixture.json").read_text(encoding="utf-8"))
TYPES = tuple(SPEC["object_types"])
LABELS = np.arange(10, 10 + len(TYPES), dtype=np.int8)


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _p39ar_context():
    frames, state, initial, thinning = _compose()
    update = run_post_thinning_standard13_membership_update(
        thinning, state["images"], SPEC["class_centroids"], state["lesion"],
        state["black"], state["outlier"], SPEC["lesion_centroid"], (0.0, 0.0),
        state["shape_values"], SPEC["shape_priors"], SPEC["smoothing_priors"],
        state["direct_relations"], frames[..., 63], state["ventricle_distances"],
        state["gray_matter_distances"], state["brainstem_distances"], SPEC["object_types"],
        centroid_options={"mode": "none", "modalities": ("T1_MPRAGE", "FLAIR"),
                          "lesion_power_exponents": (0.4, 1.2)},
        channel_weight_options={"outlier_distance": SPEC["outlier_energy"]},
        outlier_energy=SPEC["outlier_energy"],
        lesion_power_exponent=SPEC["lesion_power_exponent"], flair_channel=1,
    )
    ratio = lesion_toads_energy_ratio(initial.objective, update.membership.objective)
    input_image = nib.load(str(REGISTRATION_ROOT / "registration_inputs_ras.nii.gz"))
    registration_frames = np.asarray(input_image.dataobj, dtype=np.float32)
    shapes = tuple(registration_frames[..., 13 + index] for index in range(13))
    registered = register_after_post_thinning_membership(
        update, ratio, shapes, REGISTRATION_SPEC["initial_transform"],
        image_resolutions_mm=tuple(REGISTRATION_SPEC["image_resolutions_mm"]),
        template_resolutions_mm=tuple(REGISTRATION_SPEC["template_resolutions_mm"]),
        init_align_iter=REGISTRATION_SPEC["init_align_iter"],
    )
    assert _bits(ratio.energy_ratio) == 0xFFC00000
    assert _bits(registered.membership_update.membership.objective) == 0x7FC00000
    assert nib.aff2axcodes(input_image.affine) == ("R", "A", "S")
    return frames, update, registered, thinning.thinning


def _branch(update, registration, thinning, *, previous, iteration, n_iter_max, max_distance, source_ratio):
    return lesion_toads_post_registration_branch(
        update,
        registration.registration,
        previous_energy=previous,
        segmentation=thinning.run.segmentation,
        skeleton=thinning.run.skeleton,
        ordering=thinning.run.ordering,
        class_volume=thinning.run.class_volume,
        object_types=TYPES,
        classification=np.zeros(thinning.run.segmentation.shape, dtype=np.int16),
        template_labels=LABELS,
        lesion_class=int(SPEC["lesion_class"]),
        iteration=iteration,
        n_iter_max=n_iter_max,
        max_distance=max_distance,
        algorithm_mode="linear_fast_marching",
        relation_limit=2,
        source_energy_ratio=source_ratio,
    )


def test_p39ar_nan_and_finite_post_registration_branches_match_java() -> None:
    assert REFERENCE["source_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == REFERENCE["harness_sha256"]
    frames, update, registration, thinning = _p39ar_context()
    assert frames.shape[:3] == (17, 17, 17)
    assert len(TYPES) == 13

    cases = {case["name"]: case for case in REFERENCE["cases"]}
    actuals = {
        "nan_before_iteration_limit_negative": _branch(
            update, registration, thinning, previous=float(registration.energy_ratio.previous_energy),
            iteration=1, n_iter_max=-3, max_distance=0.1,
            source_ratio=registration.energy_ratio,
        ),
        "nan_at_iteration_limit_positive": _branch(
            update, registration, thinning, previous=float(registration.energy_ratio.previous_energy),
            iteration=3, n_iter_max=3, max_distance=0.1,
            source_ratio=registration.energy_ratio,
        ),
        "nan_at_iteration_limit_negative": _branch(
            update, registration, thinning, previous=float(registration.energy_ratio.previous_energy),
            iteration=3, n_iter_max=-3, max_distance=0.1,
            source_ratio=registration.energy_ratio,
        ),
    }
    finite_update = replace(
        update, membership=replace(update.membership, objective=np.float32(3.20))
    )
    finite_ratio = lesion_toads_energy_ratio(3.25, finite_update.membership.objective)
    finite_case_update_registration = replace(
        registration, membership_update=finite_update, energy_ratio=finite_ratio
    )
    actuals["finite_threshold_stop"] = _branch(
        finite_update, finite_case_update_registration, thinning,
        previous=3.25, iteration=1, n_iter_max=-3, max_distance=0.1,
        source_ratio=finite_ratio,
    )

    for name, result in actuals.items():
        expected = cases[name]
        assert _bits(result.energy_ratio) == int(expected["ratio_bits"], 16), name
        assert result.convergence_stop is expected["convergence_stop"], name
        assert result.iteration_limit_stop is expected["iteration_limit_stop"], name
        assert result.stop is expected["stop"], name
        assert result.run_growing is expected["run_growing"], name
        assert (result.final_classification is not None) is expected["update_classification"], name
        assert (result.relations is not None) is expected["refresh_relations"], name
        assert result.stage_order[-1] == (
            "computeRelations" if expected["refresh_relations"]
            else "updateLesionClassifcation"
        ), name
        if result.run_growing:
            assert result.growing is not None, name
            assert result.growing.events, name
            assert result.growing.remaining_queue == (), name
        else:
            assert result.growing is None, name

    with pytest.raises(ValueError, match="requires its source-computed energy ratio"):
        lesion_toads_post_registration_branch(
            update, registration.registration,
            previous_energy=float(registration.energy_ratio.previous_energy),
            segmentation=thinning.run.segmentation, skeleton=thinning.run.skeleton,
            ordering=thinning.run.ordering, class_volume=thinning.run.class_volume,
            object_types=TYPES, classification=np.zeros(thinning.run.segmentation.shape, dtype=np.int16),
            template_labels=LABELS, lesion_class=int(SPEC["lesion_class"]),
            iteration=1, n_iter_max=-3, max_distance=0.1,
            algorithm_mode="linear_fast_marching", relation_limit=2,
        )
    bad_ratio = replace(registration.energy_ratio, previous_energy=np.float32(np.inf))
    with pytest.raises(ValueError, match="source_energy_ratio must match"):
        _branch(
            update, registration, thinning,
            previous=float(registration.energy_ratio.previous_energy),
            iteration=1, n_iter_max=-3, max_distance=0.1,
            source_ratio=bad_ratio,
        )
    infinite_update = replace(
        update, membership=replace(update.membership, objective=np.float32(np.inf))
    )
    infinite_ratio = lesion_toads_energy_ratio(
        registration.energy_ratio.previous_energy, infinite_update.membership.objective
    )
    infinite_registration = replace(
        registration, membership_update=infinite_update, energy_ratio=infinite_ratio
    )
    with pytest.raises(ValueError, match="infinite source energies and ratios are unsupported"):
        _branch(
            infinite_update, infinite_registration, thinning,
            previous=float(registration.energy_ratio.previous_energy),
            iteration=1, n_iter_max=-3, max_distance=0.1,
            source_ratio=infinite_ratio,
        )
