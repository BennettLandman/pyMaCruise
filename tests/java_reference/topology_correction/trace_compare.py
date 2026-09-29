"""Compare initial labels and early Java/Python fast-marching events."""

from __future__ import annotations

from pathlib import Path
import sys
import os

import nibabel as nib
import numpy as np
from scipy import ndimage

import pymacruise.topology as topology

CASE_DIR = Path(__file__).resolve().parent
DATA_DIR = CASE_DIR.parents[1] / "data" / "synthetic" / "topology_correction"
TRACE_FILE = CASE_DIR.parents[2] / "build" / "java-reference" / "topology-correction-trace.txt"
ARGUMENTS = [argument for argument in sys.argv[1:] if argument != "--verbose"]
FIXTURE = ARGUMENTS[1] if len(ARGUMENTS) > 1 else os.environ.get("TOPOLOGY_TRACE_CASE", "spatial_hash")


def main() -> None:
    verbose = "--verbose" in sys.argv
    java_text = Path(ARGUMENTS[0]).read_text(encoding="utf-8") if ARGUMENTS else TRACE_FILE.read_text(encoding="utf-8")
    java_events = [line for line in java_text.splitlines() if line.startswith("TRACE_EVENT ")]
    java_tests = [line for line in java_text.splitlines() if line.startswith("TRACE_TEST")]
    java_decisions = [line for line in java_text.splitlines() if line.startswith("TRACE_DECISION")]
    java_total = next((line.split()[-1] for line in java_text.splitlines() if line.startswith("TRACE_EVENT_TOTAL")), "unknown")
    if verbose:
        print("JAVA INITIALIZATION")
        for line in java_text.splitlines():
            if line.startswith(("TRACE_INIT", "TRACE_LABEL", "TRACE_QUEUE_HEAD")):
                print(line)

    fixture_dir = DATA_DIR if FIXTURE == "spatial_hash" else DATA_DIR / FIXTURE
    source = nib.load(str(fixture_dir / "input.nii.gz"))
    paint = nib.load(str(fixture_dir / "paint_mask.nii.gz"))
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

    if verbose:
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
    rejected = sum(event.get("regular") is False for event in python_events)
    print(f"candidate rejected as critical: Java={sum('status=20' in line for line in java_decisions)}, Python={rejected}")

    all_events_match = len(java_events) == len(python_events)
    if all_events_match:
        for java_line, python_event in zip(java_events, python_events, strict=True):
            fields = dict(item.split("=", 1) for item in java_line.split()[2:])
            java_point = tuple(int(value) for value in fields["point"].split(","))
            java_value = float(fields["value"])
            python_value = float(python_event["value"])
            if java_point != python_event["point"] or np.float32(java_value) != np.float32(python_value):
                all_events_match = False
                break
    print(f"event trace: {'all' if all_events_match else 'NOT all'} popped points and float32 priorities match")


if __name__ == "__main__":
    main()
