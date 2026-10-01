from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_initial_loop_entry import seed_lesion_toads_initial_loop
from pymacruise.lesion_initial_thinning_composition import run_lesion_initial_thinning_pass
from pymacruise.lesion_post_thinning_membership import run_post_thinning_standard13_membership_update
from pymacruise.lesion_standard13_volume import run_lesion_toads_standard13_volume
from test_lesion_standard13_voxel import _state as voxel_state


ROOT = Path(__file__).parent / "data/synthetic/lesion_initial_thinning_composition"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))
SPEC = json.loads((ROOT / "fixture.json").read_text(encoding="utf-8"))
STATISTICS_REFERENCE = json.loads((Path(__file__).parent / "data/synthetic/lesion_post_thinning_statistics_refresh/java_reference.json").read_text(encoding="utf-8"))
UPDATE_REFERENCE = json.loads((Path(__file__).parent / "data/synthetic/lesion_post_thinning_gm_membership_update/java_reference.json").read_text(encoding="utf-8"))
HARNESS = Path(__file__).parent / "java_reference/LesionPostThinningStandard13MembershipReferenceHarness.java"


def _bits(value: np.ndarray | np.float32) -> list[int] | int:
    array = np.asarray(value, dtype=np.float32)
    if array.ndim == 0:
        return int(array.view(np.uint32))
    return array.view(np.uint32).reshape(-1).astype(np.uint64).tolist()


def _compose():
    image = nib.load(str(ROOT / "initial_thinning_inputs_ras.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    state = {
        "images": frames[..., 0:2],
        "shape_values": frames[..., 2:15],
        "direct_relations": frames[..., 15:28] > 0.5,
        "segmentation": frames[..., 28].astype(np.int16),
        "skeleton": frames[..., 29].astype(np.int16),
        "available": frames[..., 30:43] > 0.5,
        "ordering": frames[..., 43],
        "ventricle_distances": frames[..., 44].astype(np.int16),
        "gray_matter_distances": frames[..., 45].astype(np.int16),
        "brainstem_distances": frames[..., 46].astype(np.int16),
        "lesion": frames[..., 60],
        "black": frames[..., 61],
        "outlier": frames[..., 62],
    }
    initial = run_lesion_toads_standard13_volume(
        images=state["images"], class_centroids=SPEC["class_centroids"],
        class_stddev=SPEC["class_stddev"], lesion_centroid=SPEC["lesion_centroid"],
        lesion_stddev=SPEC["lesion_stddev"], channel_weights=SPEC["channel_weights"],
        shape_values=state["shape_values"], shape_priors=SPEC["shape_priors"],
        smoothing_priors=SPEC["smoothing_priors"], direct_relations=state["direct_relations"],
        initial_class_memberships=frames[..., 47:60], initial_lesion_membership=state["lesion"],
        initial_black_hole_membership=state["black"], initial_outlier_membership=state["outlier"],
        initial_lesion_relation=frames[..., 63], ventricle_distances=state["ventricle_distances"],
        gray_matter_distances=state["gray_matter_distances"], brainstem_distances=state["brainstem_distances"],
        outlier_energy=SPEC["outlier_energy"], lesion_power_exponent=SPEC["lesion_power_exponent"],
        gm_flag=False,
    )
    assert _bits(initial.class_memberships) == REFERENCE["class_membership_bits"]
    loop_entry = seed_lesion_toads_initial_loop(initial, iterations=1)
    thinning = run_lesion_initial_thinning_pass(
        loop_entry, state["segmentation"], SPEC["object_types"], state["shape_values"],
        state["skeleton"], state["available"], state["ordering"],
        lesion_class=SPEC["lesion_class"], volume_ratio=SPEC["volume_ratio"],
        limit=SPEC["thinning_limit"],
    )
    assert thinning.thinning.run.segmentation.ravel().tolist() == REFERENCE["final_segmentation"]
    return frames, state, initial, thinning


def test_post_thinning_statistics_feed_the_next_gmflag_membership_volume() -> None:
    frames, state, initial, thinning = _compose()
    result = run_post_thinning_standard13_membership_update(
        thinning,
        state["images"],
        SPEC["class_centroids"],
        state["lesion"], state["black"], state["outlier"],
        SPEC["lesion_centroid"], (0.0, 0.0),
        state["shape_values"], SPEC["shape_priors"], SPEC["smoothing_priors"],
        state["direct_relations"], frames[..., 63], state["ventricle_distances"],
        state["gray_matter_distances"], state["brainstem_distances"], SPEC["object_types"],
        centroid_options={"mode": "none", "modalities": ("T1_MPRAGE", "FLAIR"),
                          "lesion_power_exponents": (0.4, 1.2)},
        channel_weight_options={"outlier_distance": SPEC["outlier_energy"]},
        outlier_energy=SPEC["outlier_energy"],
        lesion_power_exponent=SPEC["lesion_power_exponent"], flair_channel=1,
    )
    assert result.stage_order[-4:] == (
        "computeCentroids", "computeVariances", "computeChannelWeights", "computeMemberships:GMflag=true"
    )
    stats = result.statistics
    assert _bits(stats.centroids.class_centroids) == np.asarray(STATISTICS_REFERENCE["class"], dtype=np.uint64).reshape(-1).tolist()
    assert _bits(stats.centroids.lesion_centroids) == STATISTICS_REFERENCE["aux_centroids"][0]
    assert _bits(stats.variances.class_standard_deviation) == STATISTICS_REFERENCE["class_std"]
    assert _bits(stats.variances.lesion_standard_deviation) == STATISTICS_REFERENCE["lesion_std"]
    assert _bits(stats.channel_weights) == STATISTICS_REFERENCE["weights"]

    actual = result.membership
    for key, values in {
        "class_membership_bits": actual.class_memberships,
        "class_energy_bits": actual.class_energies,
        "lesion_membership_bits": actual.lesion_membership,
        "black_hole_membership_bits": actual.black_hole_membership,
        "outlier_membership_bits": actual.outlier_membership,
        "lesion_relation_bits": actual.lesion_relation,
        "denominator_bits": actual.denominator,
        "lesion_energy_bits": actual.lesion_energy,
        "black_hole_energy_bits": actual.black_hole_energy,
        "running_objective_bits": actual.running_objective,
    }.items():
        assert _bits(values) == UPDATE_REFERENCE[key], key
    assert _bits(actual.objective) == UPDATE_REFERENCE["objective_bits"]
    assert UPDATE_REFERENCE["objective_bits"] == 0x7FC00000
    assert UPDATE_REFERENCE["class_energy_bits"].count(0x7FC00000) == 37_125
    assert UPDATE_REFERENCE["gm_flag"] is True
    assert UPDATE_REFERENCE["source_call_order"] == [
        "propagateThinning:heapLoop", "computeCentroids", "computeVariances",
        "computeChannelWeights", "computeMemberships(GMflag=true)",
    ]
    assert UPDATE_REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert UPDATE_REFERENCE["wrapper_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == UPDATE_REFERENCE["harness_sha256"]
    assert UPDATE_REFERENCE["upstream_initial_membership_reference_sha256"] == hashlib.sha256((ROOT / "java_reference.json").read_bytes()).hexdigest()
    assert hashlib.sha256((ROOT / "initial_thinning_inputs_ras.nii.gz").read_bytes()).hexdigest() == UPDATE_REFERENCE["fixture_sha256"]
    assert hashlib.sha256((Path(__file__).parent / "data/synthetic/lesion_post_thinning_statistics_refresh/java_reference.json").read_bytes()).hexdigest() == UPDATE_REFERENCE["upstream_statistics_reference_sha256"]


def test_zero_variance_source_mode_still_rejects_nonfinite_caller_inputs() -> None:
    from pymacruise.lesion_standard13_voxel import compute_lesion_toads_standard13_voxel

    state = voxel_state()
    state["image_values"] = np.asarray(state["image_values"], dtype=np.float32).copy()
    state["image_values"][0] = np.nan
    state["gm_flag"] = True
    state["source_zero_variance_nan"] = True
    with pytest.raises(ValueError, match="numeric state must contain only finite values"):
        compute_lesion_toads_standard13_voxel(**state)
