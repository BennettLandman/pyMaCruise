from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_competition import propagate_lesion_toads_competition
from pymacruise.toads_topology import toads_homeomorphic_label_change


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "lesion_competition_pass"
PASS_HARNESS = ROOT / "java_reference" / "LesionToadsCompetitionPassReferenceHarness.java"
TOPOLOGY_HARNESS = ROOT / "java_reference" / "LesionToadsCompetitionTopologyReferenceHarness.java"
JAVA_UPSTREAM = ROOT / "java_reference" / "upstream" / "edu" / "jhmi" / "rad" / "medic" / "structures"
JAVA_HEAP = JAVA_UPSTREAM / "BinaryHeap4D.java"
JAVA_LUT = JAVA_UPSTREAM / "CriticalPointLUT.java"
LUT_ASSET = ROOT.parent / "src" / "pymacruise" / "data" / "critical186LUT.raw.gz"
PASS_REFERENCE = json.loads((FIXTURE / "java_reference.json").read_text())
TOPOLOGY_REFERENCE = json.loads((FIXTURE / "java_topology_reference.json").read_text())
INPUT = FIXTURE / "competition_inputs.nii.gz"
TYPES = ("object", "object", "lesion", "blackhole")
VOLUME = np.asarray([2000.0, 1500.0, 90.0, 55.0], dtype=np.float32)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _array_hash(array: np.ndarray, dtype: str) -> str:
    little = np.asarray(array, dtype=np.dtype(dtype)).tobytes(order="C")
    return hashlib.sha256(little).hexdigest()


def _item(item) -> list[int]:
    return [
        int(item.value.view(np.uint32)), item.x, item.y, item.z, item.label,
    ]


def test_lesion_toads_topology_matches_selected_18_6_method_and_lut() -> None:
    assert TOPOLOGY_REFERENCE["source_method"] == (
        "LesionToadSegmentation.isHomeomorphicSegmentation"
    )
    assert TOPOLOGY_REFERENCE["connectivity"] == "18/6"
    assert TOPOLOGY_REFERENCE["source_sha256"] == (
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    )
    assert TOPOLOGY_REFERENCE["source_method_lines"] == "1611-1725"
    assert TOPOLOGY_REFERENCE["lut_source_sha256"] == (
        "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    )
    assert TOPOLOGY_REFERENCE["lut_asset_sha256"] == (
        "cbf1e339d78a0a2c673c939ee2c8a30ee3a35bfd9603d07ed2d48864bb0e88be"
    )
    assert _sha256(TOPOLOGY_HARNESS) == (
        "b1c820ccaa6785ee3e796c77c853a176306882456e3c02afaba38bdd5abc8d0f"
    )
    assert _sha256(JAVA_LUT) == TOPOLOGY_REFERENCE["lut_source_sha256"]
    assert _sha256(LUT_ASSET) == TOPOLOGY_REFERENCE["lut_asset_sha256"]
    assert len(TOPOLOGY_REFERENCE["cases"]) == 140
    for case in TOPOLOGY_REFERENCE["cases"]:
        neighborhood = np.fromiter(
            (int(v) for v in case["neighborhood"]), dtype=np.int16, count=27
        ).reshape((3, 3, 3))
        assert toads_homeomorphic_label_change(
            neighborhood, case["candidate"]
        ) is case["allowed"]


def test_competition_pass_matches_selected_java_on_serialized_ras_fixture() -> None:
    assert PASS_REFERENCE["source_method"] == (
        "LesionToadSegmentation.initCompetitionLabels+propagateCompetition"
    )
    assert PASS_REFERENCE["connectivity"] == "18/6"
    assert PASS_REFERENCE["source_sha256"] == (
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    )
    assert PASS_REFERENCE["source_methods"] == "674-718,2235-2327"
    assert PASS_REFERENCE["heap_sha256"] == (
        "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    )
    assert _sha256(JAVA_HEAP) == PASS_REFERENCE["heap_sha256"]
    assert _sha256(JAVA_LUT) == PASS_REFERENCE["lut_source_sha256"]
    assert _sha256(LUT_ASSET) == PASS_REFERENCE["lut_asset_sha256"]
    assert _sha256(PASS_HARNESS) == (
        "db01f237a6faa1b7a6bf05da1d636538543a41d6b3ffb82a5f5f8bb8833ceff0"
    )
    assert _sha256(INPUT) == (
        "0ad5cf2addbb7444d25cfb9bce10f6d5b1d695ddfac982602b301cbd470c6953"
    )

    image = nib.load(str(INPUT))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (0.8, 1.0, 1.2))
    packed = np.asarray(image.dataobj, dtype=np.float32)
    assert packed.shape == (17, 17, 17, 8)
    labels = np.rint(packed[..., 0]).astype(np.int16)
    memberships = packed[..., 1:5].copy()
    lesion = packed[..., 5].copy()
    blackhole = packed[..., 6].copy()
    ordering = packed[..., 7].copy()
    volume_before = VOLUME.copy()
    before = tuple(a.copy() for a in (labels, memberships, lesion, blackhole, ordering))

    for expected in PASS_REFERENCE["cases"]:
        lesion_class = expected["lesion_class"]
        actual = propagate_lesion_toads_competition(
            labels, memberships, lesion, blackhole, ordering, VOLUME,
            lesion_class=lesion_class, object_types=TYPES,
        )
        assert [_item(item) for item in actual.initialization.inserted_items] == expected["initial_queue"]
        events = [
            {
                "pop": _item(event.popped_item),
                "processed": event.processed,
                "candidate": event.candidate_label,
                "accepted": event.topology_accepted,
                "inserted": [_item(item) for item in event.inserted_neighbors],
            }
            for event in actual.events
        ]
        assert events == expected["events"]
        assert expected["remaining_queue"] == 0
        assert sum(e.processed and e.topology_accepted is True for e in actual.events) == expected["accepted"]
        assert sum(e.processed and e.topology_accepted is False for e in actual.events) == expected["rejected"]
        assert sum(not e.processed for e in actual.events) == expected["stale"]
        assert actual.remaining_queue == ()
        assert _array_hash(actual.segmentation, "<i2") == expected["final_segmentation_sha256"]
        assert _array_hash(actual.skeleton, "<i2") == expected["final_skeleton_sha256"]
        assert _array_hash(actual.available, "u1") == expected["availability_sha256"]
        assert _array_hash(actual.ordering, "<f4") == expected["ordering_sha256"]
        assert actual.class_volume.view(np.uint32).astype(np.uint64).tolist() == expected["volume_bits"]
        assert np.all(actual.initialization.available)
        np.testing.assert_array_equal(actual.initialization.skeleton, labels)
        for event in actual.events:
            item = event.popped_item
            assert all(1 <= coordinate < 16 for coordinate in (item.x, item.y, item.z))
            assert all(
                all(1 <= coordinate < 16 for coordinate in (queued.x, queued.y, queued.z))
                for queued in event.inserted_neighbors
            )
        for axis in range(3):
            for face in (0, -1):
                for actual_array, input_array in (
                    (actual.segmentation, labels),
                    (actual.skeleton, labels),
                    (actual.ordering, ordering),
                ):
                    np.testing.assert_array_equal(
                        np.take(actual_array, face, axis=axis),
                        np.take(input_array, face, axis=axis),
                    )

    for actual, expected in zip((labels, memberships, lesion, blackhole, ordering), before):
        np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(VOLUME, volume_before)
