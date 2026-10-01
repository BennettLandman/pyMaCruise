from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_standard13_voxel import compute_lesion_toads_standard13_voxel
from pymacruise.lesion_standard13_volume import run_lesion_toads_standard13_volume
from test_lesion_standard13_voxel import _state as voxel_state


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_standard13_voxel"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsStandard13VoxelReferenceHarness.java"


def _bits(value: np.ndarray | np.float32) -> list[int] | int:
    data = np.asarray(value, dtype=np.float32)
    if data.ndim == 0:
        return int(data.view(np.uint32))
    return data.view(np.uint32).reshape(-1).astype(np.uint64).tolist()


def _outputs(result: object) -> dict[str, object]:
    return {
        "class_membership_bits": _bits(result.class_memberships),
        "lesion_membership_bits": _bits(result.lesion_membership),
        "black_hole_membership_bits": _bits(result.black_hole_membership),
        "outlier_membership_bits": _bits(result.outlier_membership),
        "lesion_relation_bits": _bits(result.lesion_relation),
        "class_energy_bits": _bits(result.class_energies),
        "lesion_energy_bits": _bits(result.lesion_energy),
        "black_hole_energy_bits": _bits(result.black_hole_energy),
        "denominator_bits": _bits(result.denominator),
        "objective_bits": _bits(result.objective),
        "class_similarity_bits": _bits(result.class_similarity),
        "lesion_similarity_bits": _bits(result.lesion_similarity),
    }


def test_standard13_post_initial_gm_distance_cases_match_selected_java_harness() -> None:
    reference = json.loads((ROOT / "gm_flag_java_reference.json").read_text())
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == reference["harness_sha256"]
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    for case in reference["cases"].values():
        assert case["distance_maxima"] == {"ventricle": 2, "gray_matter": 3, "brainstem": 25}
    cases = {
        "valid_gm": (-1, 2, 10),
        "invalid_gm": (-1, -1, 10),
        "zero_gm": (-1, 0, 10),
        "zero_gm_without_brainstem": (-1, 0, -1),
        "ventricle_precedence": (1, 2, 10),
        "ventricle_only": (1, -1, 10),
        "gm_without_brainstem": (-1, 2, -1),
    }
    actual: dict[str, dict[str, object]] = {}
    for name, (vent, gm, brainstem) in cases.items():
        state = voxel_state()
        state.update(ventricle_distance=vent, gray_matter_distance=gm,
                     brainstem_distance=brainstem, gm_flag=True)
        actual[name] = _outputs(compute_lesion_toads_standard13_voxel(**state))
        for key, expected in reference["cases"][name].items():
            if key in actual[name]:
                assert actual[name][key] == expected, (name, key)

    # The Java else-if makes ventricle modulation win even when GM is valid.
    assert actual["ventricle_precedence"] == actual["ventricle_only"]
    # An invalid GM distance skips the GM branch; this equals initial mode.
    initial_mode = voxel_state()
    initial_mode.update(ventricle_distance=-1, gray_matter_distance=-1, brainstem_distance=10)
    assert actual["invalid_gm"] == _outputs(compute_lesion_toads_standard13_voxel(**initial_mode))

    # The zero-distance GM relation uses float spread/max, whereas positive GM
    # uses Java byte integer division. Brainstem is an independent later branch.
    zero_gm_relation = np.asarray(
        actual["zero_gm_without_brainstem"]["lesion_relation_bits"], dtype=np.uint32
    ).view(np.float32)
    assert zero_gm_relation == np.float32(np.float32(0.5) / np.float32(3.0))
    raw_without = np.float32(
        np.asarray(actual["gm_without_brainstem"]["lesion_membership_bits"], dtype=np.uint32).view(np.float32)
        * np.asarray(actual["gm_without_brainstem"]["denominator_bits"], dtype=np.uint32).view(np.float32)
    )
    raw_with = np.float32(
        np.asarray(actual["valid_gm"]["lesion_membership_bits"], dtype=np.uint32).view(np.float32)
        * np.asarray(actual["valid_gm"]["denominator_bits"], dtype=np.uint32).view(np.float32)
    )
    assert raw_with == np.float32(np.float32(np.float32(10.0) * raw_without) / np.float32(25.0))


def test_standard13_volume_forwards_post_initial_gm_flag_to_voxel() -> None:
    voxel = voxel_state()
    data = {"images": np.zeros((3, 3, 3, 2), dtype=np.float32)}
    data["images"][1, 1, 1] = voxel["image_values"]
    classes = np.zeros((3, 3, 3, 13), dtype=np.float32)
    classes[1, 1, 1] = voxel["class_memberships"]
    lesion = np.zeros((3, 3, 3), dtype=np.float32)
    lesion[1, 1, 1] = voxel["lesion_membership"]
    black_hole = np.zeros_like(lesion)
    outlier = np.zeros_like(lesion)
    relation = np.zeros_like(lesion)
    shape = np.broadcast_to(np.asarray(voxel["shape_values"]), classes.shape).copy()
    direct = np.broadcast_to(np.asarray(voxel["direct_relations"]), classes.shape).copy()
    faces = np.asarray(((2, 1, 1), (0, 1, 1), (1, 2, 1), (1, 0, 1), (1, 1, 2), (1, 1, 0)))
    for index, pos in enumerate(faces):
        classes[tuple(pos)] = voxel["six_face_class_memberships"][index]
        lesion[tuple(pos)] = voxel["six_face_lesion_memberships"][index]
        black_hole[tuple(pos)] = voxel["six_face_black_hole_memberships"][index]
    vent = np.full((3, 3, 3), -1, dtype=np.int16)
    gm = np.full_like(vent, -1)
    brainstem = np.full_like(vent, -1)
    vent[1, 1, 1], gm[1, 1, 1], brainstem[1, 1, 1] = -1, 0, 10
    args = {
        **data,
        "class_centroids": voxel["class_centroids"],
        "class_stddev": voxel["class_stddev"],
        "lesion_centroid": voxel["lesion_centroid"],
        "lesion_stddev": voxel["lesion_stddev"],
        "channel_weights": voxel["channel_weights"],
        "shape_values": shape,
        "shape_priors": voxel["shape_priors"],
        "smoothing_priors": voxel["smoothing_priors"],
        "direct_relations": direct,
        "initial_class_memberships": classes,
        "initial_lesion_membership": lesion,
        "initial_black_hole_membership": black_hole,
        "initial_outlier_membership": outlier,
        "initial_lesion_relation": relation,
        "ventricle_distances": vent,
        "gray_matter_distances": gm,
        "brainstem_distances": brainstem,
        "outlier_energy": voxel["outlier_energy"],
        "lesion_power_exponent": voxel["lesion_power_exponent"],
        "gm_flag": True,
    }
    volume_result = run_lesion_toads_standard13_volume(**args)
    voxel.update(ventricle_distance=-1, gray_matter_distance=0,
                 brainstem_distance=10, gm_flag=True)
    voxel_result = compute_lesion_toads_standard13_voxel(**voxel)
    for name in (
        "class_memberships", "lesion_membership", "black_hole_membership",
        "outlier_membership", "lesion_relation", "denominator", "class_energies",
        "lesion_energy", "black_hole_energy",
    ):
        np.testing.assert_array_equal(getattr(volume_result, name)[1, 1, 1], getattr(voxel_result, name))
