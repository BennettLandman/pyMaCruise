from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import lesion_toads_thinning_voxel_transition


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_thinning_transition"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
CASE_INPUTS = (
    ("direct_accept", "212221122120200120112112212", 0, False),
    ("direct_reject", "221001200120102002101100002", 0, False),
    ("fallback_accept", "222222101122202022120200001", 0, False),
    ("no_connected_target", "000000000000000000000000000", 0, False),
    ("best_already", "000000000000000000000000000", 0, False),
    ("critical_volume", "000000000000000000000000000", 0, False),
    ("lesion_volume_max", "111111111111111111111111111", 1, True),
)


def _float_bits(values) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_selected_lesion_toads_thinning_event_matches_java_excerpt() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    image = nib.load(str(ROOT / "transition_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (3, 3, 3, len(CASE_INPUTS) * 11)
    assert [case["name"] for case in REFERENCE["cases"]] == [item[0] for item in CASE_INPUTS]

    for case_index, (name, encoded, processing_class, use_max) in enumerate(CASE_INPUTS):
        base = case_index * 11
        segmentation = frames[..., base].astype(np.int16)
        expected_seg = np.fromiter((int(char) for char in encoded), dtype=np.int16).reshape(3, 3, 3)
        np.testing.assert_array_equal(segmentation, expected_seg)
        skeleton = frames[..., base + 1].astype(np.int16)
        memberships = np.asarray(
            [frames[1, 1, 1, base + 2 + c] for c in range(3)], dtype=np.float32
        )
        lesion = frames[..., base + 5]
        black_hole = frames[..., base + 6]
        shapes = np.moveaxis(frames[..., base + 7 : base + 10], -1, -1)
        volumes = np.array([12.5, 17.25, 9.0], dtype=np.float32)
        minimums = np.array([2.0, 2.0, 2.0], dtype=np.float32)
        if name == "critical_volume":
            volumes[processing_class] = minimums[processing_class]
        result = lesion_toads_thinning_voxel_transition(
            segmentation,
            skeleton,
            memberships,
            ("object", "lesion", "blackhole"),
            lesion,
            black_hole,
            shapes,
            volumes,
            minimums,
            np.array([True, True, True]),
            0.375,
            float(frames[1, 1, 1, base + 10]),
            lesion_class=1,
            processing_class=processing_class,
            use_lesion_max=use_max,
        )
        reference = REFERENCE["cases"][case_index]
        assert result.candidate_label == (None if reference["candidate"] < 0 else reference["candidate"])
        assert result.topology_accepted is reference["topology"]
        assert result.processed is reference["processed"]
        assert result.segmentation_label == reference["seg"]
        assert result.skeleton_label == reference["skeleton"]
        np.testing.assert_array_equal(_float_bits(result.class_volume), reference["volume_bits"])
        np.testing.assert_array_equal(result.available_at_center, reference["available"])
        assert int(np.asarray(result.ordering_at_center).view(np.uint32)) == (reference["ordering_bits"] & 0xFFFFFFFF)
