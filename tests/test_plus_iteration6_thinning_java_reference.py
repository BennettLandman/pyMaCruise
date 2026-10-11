from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import propagate_lesion_toads_thinning
from pymacruise.image_io import nifti_ras_on_reference, save_nifti_ras


ROOT = Path(__file__).parent / "data/synthetic/plus_iteration_6_thinning"
PROVENANCE = json.loads((ROOT / "provenance.json").read_text(encoding="utf-8"))
HARNESS = Path(__file__).parent / "java_reference/Lesion13ClassNextThinningReferenceHarness.java"
JAVA_UPSTREAM = Path(__file__).parent / "java_reference/upstream/edu/jhmi/rad/medic"
HEAP = JAVA_UPSTREAM / "structures/BinaryHeap4D.java"
LUT_JAVA = JAVA_UPSTREAM / "structures/CriticalPointLUT.java"
NUMERICS = JAVA_UPSTREAM / "utilities/Numerics.java"
LUT_GZIP = Path(__file__).parent / "../src/pymacruise/data/critical186LUT.raw.gz"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(name: str) -> np.ndarray:
    record = PROVENANCE["files"][name]
    path = ROOT / record["path"]
    assert _sha256(path) == record["sha256"]
    image = nib.load(str(path))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert list(image.shape) == record["shape"]
    return np.asarray(image.dataobj)


def _item_bytes(item) -> bytes:
    return struct.pack(">fiiii", float(item.value), item.x, item.y, item.z, item.label)


def _input_bytes(seg, skel, memberships, lesion, black, shapes, available, order) -> bytes:
    nx, ny, nz = (int(value) for value in seg.shape)
    parameters = PROVENANCE["parameters"]
    class_count = int(parameters["classes"])
    lesion_class = int(parameters["lesion_class_index"])
    type_codes = {"lesion": 1, "blackhole": 2}
    kinds = bytes(type_codes.get(value, 0) for value in parameters["object_types"])
    fields = [
        struct.pack(">5i2f", nx, ny, nz, class_count, lesion_class,
                    parameters["volume_ratio"], parameters["limit"]),
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


def _python_output_bytes(result) -> bytes:
    fields = [
        np.asarray(result.initialization.class_volume, dtype=">f4").tobytes(),
        np.asarray(result.minimum_class_volume, dtype=">f4").tobytes(),
        struct.pack(">i", len(result.initialization.inserted_items)),
    ]
    fields.extend(_item_bytes(item) for item in result.initialization.inserted_items)
    fields.append(struct.pack(">i", len(result.run.events)))
    for event in result.run.events:
        candidate = -1 if event.candidate_label is None else int(event.candidate_label)
        topology = -1 if event.topology_accepted is None else int(event.topology_accepted)
        fields.extend((
            _item_bytes(event.popped_item),
            struct.pack(">bibi", int(event.processed), candidate, topology,
                        len(event.inserted_neighbors)),
        ))
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


def _capture_python_result(
    result, reference_path: Path, output_dir: Path, remaining_queue_count: int,
) -> None:
    """Optionally save the validated state for the next source comparison."""
    reference = nib.load(str(reference_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    arrays = {
        "segmentation": result.run.segmentation.astype(np.int16, copy=False),
        "skeleton": result.run.skeleton.astype(np.int16, copy=False),
        "availability_13ch": result.run.available.astype(np.uint8, copy=False),
        "ordering": result.run.ordering.astype(np.float32, copy=False),
    }
    files = {}
    for name, values in arrays.items():
        path = output_dir / f"{name}.nii.gz"
        save_nifti_ras(nifti_ras_on_reference(reference, values), path, dtype=values.dtype)
        files[name] = {"path": path.name, "sha256": _sha256(path), "shape": list(values.shape)}
    report = {
        "case": "P41z synthetic iteration-6 source-selected thinning result",
        "validated_against": "Lesion13ClassNextThinningReferenceHarness",
        "selected_lesion_toad_source_sha256": PROVENANCE["source_call_order_sha256"]["LesionToadSegmentation.java"],
        "java_reference_harness_sha256": PROVENANCE["java_reference_harness"]["sha256"],
        "input_provenance_sha256": _sha256(ROOT / "provenance.json"),
        "synthetic_only": True,
        "spatial_convention": "RAS",
        "class_volume": [float(value) for value in result.run.class_volume],
        "minimum_class_volume": [float(value) for value in result.minimum_class_volume],
        "initial_boundary_count": len(result.initialization.inserted_items),
        "event_count": len(result.run.events),
        "remaining_queue_count": remaining_queue_count,
        "files": files,
    }
    (output_dir / "result.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def test_iteration6_thinning_matches_selected_java_reference(tmp_path: Path) -> None:
    assert PROVENANCE["input_checkpoint_iteration"] == 5
    assert PROVENANCE["synthetic_only"] is True
    assert PROVENANCE["java_reference_harness"]["sha256"] == _sha256(HARNESS)
    seg = _load("segmentation")
    skeleton = _load("skeleton")
    packed_memberships = _load("memberships_17ch")
    available = _load("availability_13ch").astype(bool)
    ordering = _load("ordering").astype(np.float32)
    shapes = _load("registered_shapes_13ch").astype(np.float32)
    assert packed_memberships.shape == (*seg.shape, 17)
    assert shapes.shape == (*seg.shape, 13)
    memberships = packed_memberships[..., :13].astype(np.float32)
    lesion = packed_memberships[..., 13].astype(np.float32)
    black_hole = packed_memberships[..., 14].astype(np.float32)
    assert all(np.isfinite(array).all() for array in (memberships, lesion, black_hole, shapes, ordering))
    parameters = PROVENANCE["parameters"]

    python_result = propagate_lesion_toads_thinning(
        seg,
        memberships,
        parameters["object_types"],
        lesion,
        black_hole,
        shapes,
        skeleton,
        available,
        ordering,
        lesion_class=parameters["lesion_class_index"],
        volume_ratio=parameters["volume_ratio"],
        limit=parameters["limit"],
        use_lesion_max=parameters["use_lesion_max"],
        # The committed inputs are explicitly checked finite above. The Java
        # NaN-specific behavior is not involved, so avoid its per-event full
        # volume scan while preserving source ordering and comparisons.
        source_faithful_nan=False,
    )
    assert python_result.initialization.inserted_items
    assert python_result.run.events
    remaining_queue_count = len(python_result.run.queue)

    javac = shutil.which("javac")
    java = shutil.which("java")
    assert javac and java, "a JDK is required for the selected Java reference harness"
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run(
        [javac, "-d", str(classes), str(HEAP), str(LUT_JAVA), str(NUMERICS), str(HARNESS)],
        check=True, capture_output=True, text=True,
    )
    resource = classes / "edu/jhmi/rad/medic/structures/critical186LUT.raw.gz"
    resource.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(LUT_GZIP, resource)
    input_path = tmp_path / "iteration6-thinning-in.bin"
    output_path = tmp_path / "iteration6-thinning-java.bin"
    input_path.write_bytes(_input_bytes(
        seg, skeleton, memberships, lesion, black_hole, shapes, available, ordering,
    ))
    subprocess.run(
        [java, "-cp", str(classes), "Lesion13ClassNextThinningReferenceHarness",
         str(input_path), str(output_path)],
        check=True, capture_output=True, text=True,
    )

    assert output_path.read_bytes() == _python_output_bytes(python_result)
    capture_dir = os.environ.get("PYMACRUISE_CAPTURE_ITERATION6_THINNING_DIR")
    if capture_dir:
        _capture_python_result(
            python_result,
            ROOT / PROVENANCE["files"]["segmentation"]["path"],
            Path(capture_dir),
            remaining_queue_count,
        )
