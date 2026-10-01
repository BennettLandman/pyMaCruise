from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
from types import SimpleNamespace

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import (
    propagate_lesion_toads_growing,
)
from pymacruise.lesion_next_iteration_thinning import run_lesion_toads_next_iteration_thinning
from pymacruise.lesion_post_registration import LesionToadsPostRegistrationBranch
from test_lesion_nan_post_registration_growth_reference import (
    HARNESS as GROW_HARNESS,
    HEAP,
    LUT_GZIP,
    LUT_JAVA,
    NUMERICS,
    _expected_bytes as expected_growth_bytes,
    _input_bytes as growth_input_bytes,
)


ROOT = Path(__file__).parent / "data/synthetic/lesion13_finite_grow_thin"
SPEC = json.loads((ROOT / "fixture.json").read_text(encoding="utf-8"))
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))
THIN_HARNESS = Path(__file__).parent / "java_reference/Lesion13ClassNextThinningReferenceHarness.java"
TYPES = tuple(SPEC["object_types"])
CLASSES = int(SPEC["class_count"])
LESION = int(SPEC["lesion_class"])


def _item_bytes(item) -> bytes:
    return struct.pack(">fiiii", float(item.value), item.x, item.y, item.z, item.label)


def _thinning_input_bytes(seg, skel, memberships, lesion, black, shapes, available, order) -> bytes:
    nx, ny, nz = (int(value) for value in seg.shape)
    kinds = bytes(1 if value == "lesion" else 2 if value == "blackhole" else 0 for value in TYPES)
    fields = [
        struct.pack(">5i2f", nx, ny, nz, CLASSES, LESION, 0.4, 1000.0),
        kinds,
        np.asarray(seg, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(skel, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(memberships, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(lesion, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(black, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(shapes, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(available, dtype=np.uint8).ravel(order="C").tobytes(),
        np.asarray(order, dtype=">f4").ravel(order="C").tobytes(),
    ]
    return b"".join(fields)


def _expected_thinning_bytes(result) -> bytes:
    fields = [
        np.asarray(result.initialization.class_volume, dtype=">f4").tobytes(),
        np.asarray(result.minimum_class_volume, dtype=">f4").tobytes(),
        struct.pack(">i", len(result.initialization.inserted_items)),
    ]
    fields.extend(_item_bytes(item) for item in result.initialization.inserted_items)
    fields.append(struct.pack(">i", len(result.run.events)))
    for event in result.run.events:
        fields.extend((_item_bytes(event.popped_item), struct.pack(">bi", int(event.processed), len(event.inserted_neighbors))))
        fields.extend(_item_bytes(item) for item in event.inserted_neighbors)
    fields.extend((
        np.asarray(result.run.segmentation, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(result.run.skeleton, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(result.run.available, dtype=np.uint8).ravel(order="C").tobytes(),
        np.asarray(result.run.ordering, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(result.run.class_volume, dtype=">f4").tobytes(),
    ))
    remaining = []
    while result.run.queue:
        remaining.append(result.run.queue.pop_first())
    fields.append(struct.pack(">i", len(remaining)))
    fields.extend(_item_bytes(item) for item in remaining)
    return b"".join(fields)


def test_finite_13class_state_growth_then_next_thinning_matches_java(tmp_path: Path) -> None:
    assert REFERENCE["lesion_toads_source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert REFERENCE["wrapper_source_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    assert hashlib.sha256(GROW_HARNESS.read_bytes()).hexdigest() == REFERENCE["grow_harness_sha256"]
    assert hashlib.sha256(THIN_HARNESS.read_bytes()).hexdigest() == REFERENCE["thin_harness_sha256"]
    assert hashlib.sha256(HEAP.read_bytes()).hexdigest() == REFERENCE["binary_heap_sha256"]
    assert hashlib.sha256(NUMERICS.read_bytes()).hexdigest() == REFERENCE["numerics_sha256"]
    image_path = ROOT / "finite_13class_grow_thin_ras.nii.gz"
    assert hashlib.sha256(image_path.read_bytes()).hexdigest() == SPEC["sha256_nii_gz"]
    image = nib.load(str(image_path))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    packed = np.asarray(image.dataobj, dtype=np.float32)
    assert hashlib.sha256(packed.astype("<f4", copy=False).tobytes()).hexdigest() == SPEC["sha256_uncompressed_float32le"]
    seg = packed[..., 0].astype(np.int16)
    skel = packed[..., 1].astype(np.int16)
    memberships = packed[..., 2:15]
    lesion, black = packed[..., 15], packed[..., 16]
    shapes, order = packed[..., 17:30], packed[..., 30]
    assert memberships.shape[-1] == len(TYPES) == CLASSES == 13
    assert np.isfinite(packed).all()
    assert tuple(packed.shape[:3]) == tuple(SPEC["shape_xyz"])
    support = np.bincount(seg[1:-1, 1:-1, 1:-1].ravel(), minlength=CLASSES)
    assert np.all(support > 0)
    center = tuple(SPEC["growth_center_xyz"])
    assert skel[center] == -1

    class_volume = np.zeros(CLASSES, dtype=np.float32)
    interior_seg = seg[1:-1, 1:-1, 1:-1]
    interior_memberships = memberships[1:-1, 1:-1, 1:-1]
    for k in range(CLASSES):
        values = interior_memberships[..., k][interior_seg == k]
        if k == LESION:
            values = values + lesion[1:-1, 1:-1, 1:-1][interior_seg == k]
            values = values + black[1:-1, 1:-1, 1:-1][interior_seg == k]
        class_volume[k] = np.sum(values, dtype=np.float32)
    grown = propagate_lesion_toads_growing(
        seg, skel, memberships, order, class_volume, TYPES, lesion, black,
        lesion_class=LESION,
    )
    assert len(grown.events) == 3
    assert sum(event.topology_accepted is True for event in grown.events) == 1
    assert len(grown.events) == REFERENCE["growth_event_count"]
    assert sum(event.topology_accepted is True for event in grown.events) == REFERENCE["growth_accepted_count"]
    assert grown.segmentation[center] == SPEC["growth_center_class"]
    assert grown.remaining_queue == ()

    classes_dir = tmp_path / "classes"
    classes_dir.mkdir()
    javac = shutil.which("javac")
    java = shutil.which("java")
    assert javac and java, "a JDK is required for the pinned Java method references"
    subprocess.run(
        [javac, "-d", str(classes_dir), str(HEAP), str(LUT_JAVA), str(NUMERICS),
         str(GROW_HARNESS), str(THIN_HARNESS)],
        check=True, capture_output=True, text=True,
    )
    resource = classes_dir / "edu/jhmi/rad/medic/structures/critical186LUT.raw.gz"
    resource.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(LUT_GZIP, resource)
    growth_input, growth_output = tmp_path / "grow-in.bin", tmp_path / "grow-out.bin"
    growth_input.write_bytes(growth_input_bytes(seg, skel, memberships, order, class_volume, lesion, black))
    subprocess.run([java, "-cp", str(classes_dir), "LesionGrowingPassDynamicReferenceHarness", str(growth_input), str(growth_output)], check=True, capture_output=True, text=True)
    assert growth_output.read_bytes() == expected_growth_bytes(grown)

    continued_branch = LesionToadsPostRegistrationBranch(
        membership_result=SimpleNamespace(
            membership=SimpleNamespace(
                class_memberships=memberships,
                lesion_membership=lesion,
                black_hole_membership=black,
            )
        ),
        registration_handoff=SimpleNamespace(),
        energy_ratio=np.float32(0.5),
        convergence_stop=False,
        iteration_limit_stop=False,
        stop=False,
        stop_reason=None,
        run_growing=True,
        growing=grown,
        segmentation=grown.segmentation,
        final_classification=None,
        relations=SimpleNamespace(),
        stage_order=("propagateGrowing", "computeRelations"),
    )
    thinned = run_lesion_toads_next_iteration_thinning(
        continued_branch,
        shapes,
        object_types=TYPES,
        lesion_class=LESION,
        volume_ratio=0.4,
        limit=1000.0,
    )
    assert len(thinned.initialization.inserted_items) == REFERENCE["thinning_initial_queue_count"]
    assert len(thinned.run.events) == REFERENCE["thinning_event_count"]
    assert sum(event.processed for event in thinned.run.events) == REFERENCE["thinning_processed_count"]
    thin_input, thin_output = tmp_path / "thin-in.bin", tmp_path / "thin-out.bin"
    thin_input.write_bytes(_thinning_input_bytes(
        continued_branch.segmentation, grown.skeleton, memberships, lesion, black, shapes,
        grown.available, grown.ordering,
    ))
    subprocess.run([java, "-cp", str(classes_dir), "Lesion13ClassNextThinningReferenceHarness", str(thin_input), str(thin_output)], check=True, capture_output=True, text=True)
    assert thin_output.read_bytes() == _expected_thinning_bytes(thinned)
    assert len(thinned.run.events) > 0
