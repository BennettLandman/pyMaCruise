from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
import struct
import subprocess

import numpy as np

from test_lesion_nan_post_registration_branch import _branch, _p39ar_context
from test_lesion_post_thinning_gm_membership import SPEC


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference/LesionGrowingPassDynamicReferenceHarness.java"
HEAP = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java"
LUT_JAVA = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/structures/CriticalPointLUT.java"
NUMERICS = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/utilities/Numerics.java"
LUT_GZIP = ROOT.parent / "src/pymacruise/data/critical186LUT.raw.gz"
HARNESS_SHA256 = "e324c89358dc83096c2cf98ef9a64f404263c0566292f25021b2b74bdfab7d62"
TYPES = tuple(SPEC["object_types"])


def _input_bytes(segmentation, skeleton, memberships, ordering, volumes, lesion, black_hole) -> bytes:
    shape = tuple(int(v) for v in segmentation.shape)
    classes = int(memberships.shape[-1])
    kinds = bytes(1 if name == "lesion" else 2 if name == "blackhole" else 0 for name in TYPES)
    fields = [
        struct.pack(">6i", *shape, classes, int(SPEC["lesion_class"]), 0),
        kinds,
        np.asarray(segmentation, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(skeleton, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(memberships, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(ordering, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(volumes, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(lesion, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(black_hole, dtype=">f4").ravel(order="C").tobytes(),
    ]
    return b"".join(fields)


def _item_bytes(item) -> bytes:
    return struct.pack(">fiiii", float(item.value), item.x, item.y, item.z, item.label)


def _expected_bytes(growing) -> bytes:
    fields = [
        np.asarray(growing.segmentation, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(growing.skeleton, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(growing.available, dtype=np.uint8).ravel(order="C").tobytes(),
        np.asarray(growing.ordering, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(growing.class_volume, dtype=">f4").ravel(order="C").tobytes(),
        struct.pack(">i", len(growing.events)),
    ]
    for event in growing.events:
        accepted = -1 if event.topology_accepted is None else int(event.topology_accepted)
        fields.extend((_item_bytes(event.popped_item), struct.pack(">bI", accepted, len(event.inserted_neighbors))))
        fields.extend(_item_bytes(item) for item in event.inserted_neighbors)
    return b"".join(fields)


def test_p39ar_growth_on_folded_13_class_state_matches_java(tmp_path: Path) -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    javac, java = shutil.which("javac"), shutil.which("java")
    assert javac is not None and java is not None, "a JDK is required for pinned Java references"

    frames, update, registration, thinning = _p39ar_context()
    before = thinning.run
    result = _branch(
        update, registration, thinning,
        previous=float(registration.energy_ratio.previous_energy), iteration=3,
        n_iter_max=3, max_distance=0.1, source_ratio=registration.energy_ratio,
    )
    assert result.run_growing and result.growing is not None
    assert frames.shape[:3] == before.segmentation.shape == (17, 17, 17)
    assert len(TYPES) == 13
    assert len(result.growing.events) == 1860

    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HEAP), str(LUT_JAVA), str(NUMERICS), str(HARNESS)], check=True, capture_output=True, text=True)
    resource = classes / "edu/jhmi/rad/medic/structures/critical186LUT.raw.gz"
    resource.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(LUT_GZIP, resource)

    source = tmp_path / "growth-input.bin"
    output = tmp_path / "growth-output.bin"
    source.write_bytes(_input_bytes(
        before.segmentation, before.skeleton,
        result.membership_result.membership.class_memberships,
        before.ordering, before.class_volume,
        result.membership_result.membership.lesion_membership,
        result.membership_result.membership.black_hole_membership,
    ))
    subprocess.run([java, "-cp", str(classes), "LesionGrowingPassDynamicReferenceHarness", str(source), str(output)], check=True, capture_output=True, text=True)
    actual = output.read_bytes()
    expected = _expected_bytes(result.growing)
    assert actual == expected
