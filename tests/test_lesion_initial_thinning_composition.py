from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_initial_loop_entry import seed_lesion_toads_initial_loop
from pymacruise.lesion_initial_thinning_composition import run_lesion_initial_thinning_pass
from pymacruise.lesion_standard13_volume import run_lesion_toads_standard13_volume


ROOT = Path(__file__).parent / "data/synthetic/lesion_initial_thinning_composition"
HARNESS = Path(__file__).parent / "java_reference/LesionToadsInitialThinningCompositionReferenceHarness.java"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))
SPEC = json.loads((ROOT / "fixture.json").read_text(encoding="utf-8"))
CLASSES = 13


def _unsigned_bits(values: np.ndarray | np.float32 | float) -> list[int] | int:
    array = np.asarray(values, dtype=np.float32)
    bits = array.view(np.uint32)
    if bits.ndim == 0:
        return int(bits)
    return bits.reshape(-1).astype(np.uint64).tolist()


def _signed_bits(value: np.float32 | float) -> int:
    bits = int(np.asarray(value, dtype=np.float32).view(np.uint32))
    return bits if bits < 0x80000000 else bits - 0x100000000


def _source_state() -> dict[str, object]:
    image = nib.load(str(ROOT / "initial_thinning_inputs_ras.nii.gz"))
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert image.shape == tuple(SPEC["shape"])
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert np.array_equal(image.affine, np.eye(4))
    assert _unsigned_bits(frames[..., 0:2]) is not None

    return {
        "images": frames[..., 0:2],
        "shape_values": frames[..., 2:15],
        "direct_relations": frames[..., 15:28] > 0.5,
        "segmentation": frames[..., 28].astype(np.int16),
        "initial_skeleton": frames[..., 29].astype(np.int16),
        "initial_available": frames[..., 30:43] > 0.5,
        "initial_ordering": frames[..., 43],
        "ventricle_distances": frames[..., 44].astype(np.int16),
        "gray_matter_distances": frames[..., 45].astype(np.int16),
        "brainstem_distances": frames[..., 46].astype(np.int16),
        "initial_class_memberships": frames[..., 47:60],
        "initial_lesion_membership": frames[..., 60],
        "initial_black_hole_membership": frames[..., 61],
        "initial_outlier_membership": frames[..., 62],
        "initial_lesion_relation": frames[..., 63],
    }


def test_initial_gm_false_memberships_feed_first_selected_thinning_pass() -> None:
    image = nib.load(str(ROOT / "initial_thinning_inputs_ras.nii.gz"))
    packed = np.asarray(image.dataobj, dtype=np.float32)
    packed_hash = hashlib.sha256(packed.astype("<f4", copy=False).tobytes()).hexdigest()
    image_hash = hashlib.sha256(packed[..., :2].astype("<f4", copy=False).tobytes()).hexdigest()
    assert packed_hash == SPEC["packed_sha256_float32le"]
    assert image_hash == SPEC["image_sha256_float32le"] == REFERENCE["image_sha256"]
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == "32985612ee9ef1487f7fc170206f5cd11cb826f5f9873c0955dfec4380fc2d57"
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert REFERENCE["wrapper_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    assert REFERENCE["power_table_source_sha256"] == "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
    assert REFERENCE["power_table_float_exponent_bits"] == "3f99999a"
    assert REFERENCE["heap_source_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert REFERENCE["topology_lut_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert REFERENCE["source_call_order"] == [
        "computeMemberships(GMflag=false)",
        "setGMflag",
        "Eprev=initial_energy",
        "Eratio=1.0f",
        "initThinningLabels",
        "propagateThinning",
    ]

    state = _source_state()
    names = (
        "class_centroids", "class_stddev", "lesion_centroid", "lesion_stddev",
        "channel_weights", "shape_priors", "smoothing_priors",
    )
    run_state = {
        **{name: np.asarray(SPEC[name], dtype=np.float32) for name in names},
        **{name: state[name] for name in (
            "images", "shape_values", "direct_relations", "initial_class_memberships",
            "initial_lesion_membership", "initial_black_hole_membership",
            "initial_outlier_membership", "initial_lesion_relation", "ventricle_distances",
            "gray_matter_distances", "brainstem_distances",
        )},
        "outlier_energy": SPEC["outlier_energy"],
        "lesion_power_exponent": SPEC["lesion_power_exponent"],
        "gm_flag": False,
    }
    fixed_inputs_before = {
        name: np.asarray(value).copy()
        for name, value in run_state.items()
        if isinstance(value, np.ndarray)
    }
    initial = run_lesion_toads_standard13_volume(**run_state)
    entry = seed_lesion_toads_initial_loop(initial, iterations=1)
    assert entry.gm_flag is True
    assert entry.initial_membership is initial
    assert entry.next_stage == "propagateThinning"

    result = run_lesion_initial_thinning_pass(
        entry,
        state["segmentation"],
        SPEC["object_types"],
        state["shape_values"],
        state["initial_skeleton"],
        state["initial_available"],
        state["initial_ordering"],
        lesion_class=SPEC["lesion_class"],
        volume_ratio=SPEC["volume_ratio"],
        limit=SPEC["thinning_limit"],
        use_lesion_max=SPEC["use_lesion_max"],
    )
    assert result.stage_order[-2:] == (
        "propagateThinning:initThinningLabels",
        "propagateThinning:heapLoop",
    )
    assert not result.thinning.run.queue

    initial_maps = {
        "class_membership_bits": initial.class_memberships,
        "class_energy_bits": initial.class_energies,
        "lesion_membership_bits": initial.lesion_membership,
        "black_hole_membership_bits": initial.black_hole_membership,
        "outlier_membership_bits": initial.outlier_membership,
        "lesion_relation_bits": initial.lesion_relation,
        "denominator_bits": initial.denominator,
        "lesion_energy_bits": initial.lesion_energy,
        "black_hole_energy_bits": initial.black_hole_energy,
        "running_objective_bits": initial.running_objective,
    }
    for name, values in initial_maps.items():
        assert _unsigned_bits(values) == REFERENCE[name], name
    assert _unsigned_bits(initial.objective) == REFERENCE["objective_bits"]

    frontier = result.thinning.initialization
    assert _unsigned_bits(frontier.class_volume) == REFERENCE["initial_volume_bits"]
    assert _unsigned_bits(result.thinning.minimum_class_volume) == REFERENCE["minimum_volume_bits"]
    assert _unsigned_bits(result.thinning.run.class_volume) == REFERENCE["final_volume_bits"]

    inserted = [
        [_signed_bits(item.value), item.x, item.y, item.z, item.label]
        for item in frontier.inserted_items
    ]
    assert inserted == REFERENCE["initial_queue"]
    events = []
    for event in result.thinning.run.events:
        neighbors = [
            [_signed_bits(item.value), item.x, item.y, item.z, item.label]
            for item in event.inserted_neighbors
        ]
        events.append([
            _signed_bits(event.popped_item.value), event.popped_item.x,
            event.popped_item.y, event.popped_item.z, event.popped_item.label,
            event.processed, event.candidate_label, event.topology_accepted, neighbors,
        ])
    assert events == REFERENCE["events"]
    assert len(events) == 1917
    assert sum(event[5] for event in events) == 1247
    assert sum(event[7] is True for event in events) == 626
    assert sum(event[7] is False for event in events) == 299

    np.testing.assert_array_equal(result.thinning.run.segmentation.ravel(), REFERENCE["final_segmentation"])
    np.testing.assert_array_equal(result.thinning.run.skeleton.ravel(), REFERENCE["final_skeleton"])
    np.testing.assert_array_equal(result.thinning.run.available.ravel(), REFERENCE["final_available"])
    assert _unsigned_bits(result.thinning.run.ordering) == REFERENCE["final_ordering_bits"]
    for name, before in fixed_inputs_before.items():
        np.testing.assert_array_equal(run_state[name], before, err_msg=f"mutated fixed input {name}")


def test_first_thinning_rejects_non_loop_or_post_membership_state() -> None:
    packed = np.asarray(
        nib.load(str(ROOT / "initial_thinning_inputs_ras.nii.gz")).dataobj,
        dtype=np.float32,
    )
    state = _source_state()
    initials = run_lesion_toads_standard13_volume(
        images=state["images"],
        class_centroids=SPEC["class_centroids"],
        class_stddev=SPEC["class_stddev"],
        lesion_centroid=SPEC["lesion_centroid"],
        lesion_stddev=SPEC["lesion_stddev"],
        channel_weights=SPEC["channel_weights"],
        shape_values=state["shape_values"],
        shape_priors=SPEC["shape_priors"],
        smoothing_priors=SPEC["smoothing_priors"],
        direct_relations=state["direct_relations"],
        initial_class_memberships=state["initial_class_memberships"],
        initial_lesion_membership=state["initial_lesion_membership"],
        initial_black_hole_membership=state["initial_black_hole_membership"],
        initial_outlier_membership=state["initial_outlier_membership"],
        initial_lesion_relation=state["initial_lesion_relation"],
        ventricle_distances=state["ventricle_distances"],
        gray_matter_distances=state["gray_matter_distances"],
        brainstem_distances=state["brainstem_distances"],
        outlier_energy=SPEC["outlier_energy"],
        lesion_power_exponent=SPEC["lesion_power_exponent"],
        gm_flag=False,
    )
    stopped = seed_lesion_toads_initial_loop(initials, iterations=0)
    try:
        run_lesion_initial_thinning_pass(
            stopped, state["segmentation"], SPEC["object_types"], state["shape_values"],
            state["initial_skeleton"], state["initial_available"], state["initial_ordering"],
            lesion_class=SPEC["lesion_class"], volume_ratio=SPEC["volume_ratio"],
            limit=SPEC["thinning_limit"],
        )
    except ValueError as exc:
        assert "does not reach propagateThinning" in str(exc)
    else:
        raise AssertionError("non-loop state unexpectedly reached thinning")
    assert packed.shape == tuple(SPEC["shape"])
