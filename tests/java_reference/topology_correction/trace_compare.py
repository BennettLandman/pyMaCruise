"""Compare initial labels and early Java/Python fast-marching events."""

from __future__ import annotations

from pathlib import Path
import sys

import nibabel as nib
import numpy as np
from scipy import ndimage

import pymacruise.topology as topology

CASE_DIR = Path(__file__).resolve().parent
DATA_DIR = CASE_DIR.parents[1] / "data" / "synthetic" / "topology_correction"
TRACE_FILE = CASE_DIR.parents[2] / "build" / "java-reference" / "topology-correction-trace.txt"


def main() -> None:
    java_text = Path(sys.argv[1]).read_text(encoding="utf-8") if len(sys.argv) > 1 else TRACE_FILE.read_text(encoding="utf-8")
    java_events = [line for line in java_text.splitlines() if line.startswith("TRACE_EVENT ")]
    java_tests = [line for line in java_text.splitlines() if line.startswith("TRACE_TEST")]
    java_decisions = [line for line in java_text.splitlines() if line.startswith("TRACE_DECISION")]
    java_total = next((line.split()[-1] for line in java_text.splitlines() if line.startswith("TRACE_EVENT_TOTAL")), "unknown")
    print("JAVA INITIALIZATION")
    for line in java_text.splitlines():
        if line.startswith(("TRACE_INIT", "TRACE_LABEL", "TRACE_QUEUE_HEAD")):
            print(line)

    source = nib.load(str(DATA_DIR / "input.nii.gz"))
    paint = nib.load(str(DATA_DIR / "paint_mask.nii.gz"))
    python_events: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    original_lower = topology._lowest_object_neighbor
    original_regular = topology._is_regular_object_addition

    def trace_lower(distance, labels, point, value):
        nonlocal current
        lower = original_lower(distance, labels, point, value)
        current = {
            "point": tuple(int(v) for v in point),
            "value": float(value),
            "lower": float(lower),
            "label": int(labels[point]),
        }
        python_events.append(current)
        return lower

    def trace_regular(labels, point):
        neighborhood = labels[
            point[0] - 1 : point[0] + 2,
            point[1] - 1 : point[1] + 2,
            point[2] - 1 : point[2] + 2,
        ]
        candidate = neighborhood >= topology._TEST
        foreground, background = topology.topology_numbers_618(candidate)
        source_labels, _ = ndimage.label(candidate, structure=topology._FOREGROUND_6)
        source_fg = len(
            {
                int(source_labels[face])
                for face in ((0, 1, 1), (2, 1, 1), (1, 0, 1), (1, 2, 1), (1, 1, 0), (1, 1, 2))
                if source_labels[face] != 0
            }
        )
        regular = original_regular(labels, point)
        if current is not None:
            current["foreground"] = foreground
            current["source_foreground"] = source_fg
            current["background"] = background
            current["regular"] = regular
        return regular

    topology._lowest_object_neighbor = trace_lower
    topology._is_regular_object_addition = trace_regular
    try:
        topology.topology_correct(
            source,
            paint,
            paint_threshold=0.5,
            lowest_level=0.0,
            highest_level=1.0,
            regularize=0.05,
        )
    finally:
        topology._lowest_object_neighbor = original_lower
        topology._is_regular_object_addition = original_regular

    print("\nFIRST 24 QUEUE EVENTS")
    print("Java event format: event point value lower; then foreground/background/regular decision when tested")
    for index in range(min(24, len(python_events))):
        event = python_events[index]
        print(
            f"PY_EVENT {index} point={event['point']} value={event['value']:.9g} "
            f"lower={event['lower']:.9g} label={event['label']} "
            f"fg={event.get('foreground','-')} java-style-fg={event.get('source_foreground','-')} "
            f"bg={event.get('background','-')} "
            f"regular={event.get('regular','not-tested')}"
        )
        if index < len(java_events):
            print("  " + java_events[index])
            matching = next((line for line in java_tests if line.split()[1] == str(index)), None)
            if matching:
                print("  " + matching)
            decision = next((line for line in java_decisions if line.split()[1] == str(index)), None)
            if decision:
                print("  " + decision)
    print(f"complete queue pops: Java={java_total}, Python={len(python_events)}")

    java_event_records: list[tuple[tuple[int, int, int], np.float32]] = []
    for line in java_events:
        fields = line.split()
        point = tuple(int(value) for value in fields[2].split("=")[1].split(","))
        value = np.float32(fields[3].split("=")[1])
        java_event_records.append((point, value))
    first_divergence = None
    for index, (java_record, python_record) in enumerate(zip(java_event_records, python_events)):
        java_point, java_value = java_record
        python_point = python_record["point"]
        python_value = np.float32(python_record["value"])
        if java_point != python_point:
            kind = "priority/tie ordering" if java_value == python_value else "queue ordering after different float32 priorities"
            first_divergence = (index, kind, java_point, float(java_value), python_point, float(python_value))
            break
        if java_value != python_value:
            first_divergence = (index, "float32 candidate arithmetic", java_point, float(java_value), python_point, float(python_value))
            break
    if first_divergence:
        index, kind, java_point, java_value, python_point, python_value = first_divergence
        print(
            f"first divergence: event {index}, {kind}; Java={java_point} @ {java_value:.9g}, "
            f"Python={python_point} @ {python_value:.9g}"
        )
    elif len(java_event_records) != len(python_events):
        print(f"first divergence: queue exhaustion length Java={len(java_event_records)}, Python={len(python_events)}")
    else:
        print("event trace: all popped points and float32 priorities match")


if __name__ == "__main__":
    main()
