from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import propagate_lesion_toads_growing


ROOT = Path(__file__).parent / "data/synthetic/lesion_growing_pass"
REFERENCE = ROOT / "java_reference.json"
N = 9
CLASSES = 4
FACE = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _skeleton() -> np.ndarray:
    skeleton = np.full((N, N, N), -1, dtype=np.int16)
    skeleton[3, 4, 3] = 1
    skeleton[5, 2, 5] = 2
    skeleton[1, 4, 7] = 3
    return skeleton


def _memberships() -> tuple[np.ndarray, ...]:
    memberships = np.empty((N, N, N, CLASSES), dtype=np.float32)
    lesion = np.empty((N, N, N), dtype=np.float32)
    black_hole = np.empty((N, N, N), dtype=np.float32)
    ordering = np.empty((N, N, N), dtype=np.float32)
    for x in range(N):
        for y in range(N):
            for z in range(N):
                ordering[x, y, z] = -((x * 13 + y * 7 + z * 3) % 19) / 10.0
                lesion[x, y, z] = ((x + 3 * y + 5 * z) % 7 + 1) / 10.0
                black_hole[x, y, z] = ((2 * x + y + 3 * z) % 9 + 1) / 12.0
                for k in range(CLASSES):
                    memberships[x, y, z, k] = ((x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1) / 30.0
    return memberships, lesion, black_hole, ordering


def _labels(text: str, *, empty: bool = False) -> np.ndarray:
    decode = lambda item: -1 if empty and item == "a" else int(item)
    return np.fromiter((decode(item) for item in text), dtype=np.int16).reshape((N, N, N))


def _event_trace(events) -> list[list[object]]:
    trace = []
    for event in events:
        item = event.popped_item
        popped = [_bits(item.value), item.x, item.y, item.z, item.label]
        topology = -1 if event.topology_accepted is None else int(event.topology_accepted)
        inserted = [
            [_bits(neighbor.value), neighbor.x, neighbor.y, neighbor.z, neighbor.label]
            for neighbor in event.inserted_neighbors
        ]
        trace.append([popped, topology, inserted])
    return trace


def test_folded_lesion_growing_pass_matches_java_through_queue_exhaustion() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert reference["lut_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"

    image = nib.load(str(ROOT / "folded_labels.nii.gz"))
    assert image.shape == (N, N, N)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    segmentation = np.asarray(image.dataobj, dtype=np.int16)
    skeleton = _skeleton()
    memberships, lesion, black_hole, ordering = _memberships()
    assert np.array_equal(segmentation, _labels(reference["cases"][0]["initial_segmentation"]))

    assert [case["name"] for case in reference["cases"]] == ["sum", "max"]
    for case in reference["cases"]:
        assert np.array_equal(skeleton, _labels(case["initial_skeleton"], empty=True))
        result = propagate_lesion_toads_growing(
            segmentation,
            skeleton,
            memberships,
            ordering,
            np.asarray([1, 2, 3, 4], dtype=np.float32),
            ("obj", "lesion", "blackhole", "optimized_lesion"),
            lesion,
            black_hole,
            lesion_class=1,
            use_lesion_max=bool(case["use_lesion_max"]),
        )

        assert _event_trace(result.events) == case["events"], case["name"]
        assert np.array_equal(result.segmentation, _labels(case["segmentation"]))
        assert np.array_equal(result.skeleton, _labels(case["skeleton"], empty=True))
        assert "".join("1" if value else "0" for value in result.available.ravel()) == case["available"]
        assert [_bits(value) for value in result.ordering.ravel()] == case["ordering_bits"]
        assert [_bits(value) for value in result.class_volume] == case["volume_bits"]
        assert result.remaining_queue == ()
        assert result.events
        assert any(event.topology_accepted is False for event in result.events)
        assert any(event.topology_accepted is None for event in result.events)

        compact = propagate_lesion_toads_growing(
            segmentation,
            skeleton,
            memberships,
            ordering,
            np.asarray([1, 2, 3, 4], dtype=np.float32),
            ("obj", "lesion", "blackhole", "optimized_lesion"),
            lesion,
            black_hole,
            lesion_class=1,
            use_lesion_max=bool(case["use_lesion_max"]),
            retain_event_states=False,
        )
        assert _event_trace(compact.events) == case["events"], case["name"]
        assert np.array_equal(compact.segmentation, result.segmentation)
        assert np.array_equal(compact.skeleton, result.skeleton)
        assert np.array_equal(compact.available, result.available)
        assert np.array_equal(compact.ordering, result.ordering)
        assert np.array_equal(compact.class_volume, result.class_volume)
        assert compact.events
        assert all(not hasattr(event, "segmentation") for event in compact.events)

        # MaCRUISE's standard-13 descriptor stores uppercase names. They must
        # select the same lesion scoring as the lowercase source type names.
        uppercase = propagate_lesion_toads_growing(
            segmentation,
            skeleton,
            memberships,
            ordering,
            np.asarray([1, 2, 3, 4], dtype=np.float32),
            tuple(value.upper() for value in ("obj", "lesion", "blackhole", "optimized_lesion")),
            lesion,
            black_hole,
            lesion_class=1,
            use_lesion_max=bool(case["use_lesion_max"]),
            retain_event_states=False,
        )
        assert _event_trace(uppercase.events) == _event_trace(compact.events)
        assert np.array_equal(uppercase.segmentation, compact.segmentation)
        assert np.array_equal(uppercase.skeleton, compact.skeleton)
        assert np.array_equal(uppercase.available, compact.available)
        assert np.array_equal(uppercase.ordering, compact.ordering)
        assert np.array_equal(uppercase.class_volume, compact.class_volume)

    assert not np.array_equal(reference["cases"][0]["volume_bits"], reference["cases"][1]["volume_bits"])
