#!/usr/bin/env python3
"""Replay the selected Java/Python iteration-5→6 core seams in order.

The committed iteration-5 checkpoint supplies the thinning inputs. Its
source-selected thinning result is then tied byte-for-byte to the committed
post-thinning state consumed by the selected statistics and membership
comparators. This is a composition of selected Java method excerpts, not an
independent complete Java algorithm run.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(TESTS))
sys.path.insert(0, str(ROOT / "scripts"))

import test_plus_iteration6_thinning_java_reference as thinning_test  # noqa: E402
from pymacruise.lesion_growth import propagate_lesion_toads_thinning  # noqa: E402
from compare_p41z_next_iteration_statistics_java import compare as compare_statistics  # noqa: E402
from compare_p41z_iteration6_membership_java import compare as compare_membership  # noqa: E402

MEMBERSHIP_CASE = TESTS / "data/synthetic/p41z_iteration6_membership_java"
THINNING_RESULT = TESTS / "data/synthetic/plus_iteration_6_thinning_result"
RESULT_PATH = ROOT / "docs/reports/p41z-iteration5-to6-core-replay-2026-10-11.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_ras(path: Path) -> np.ndarray:
    image = nib.load(str(path))
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ValueError(f"not RAS: {path}")
    return np.asarray(image.dataobj)


def _run_thinning() -> dict:
    provenance = thinning_test.PROVENANCE
    seg = thinning_test._load("segmentation")
    skeleton = thinning_test._load("skeleton")
    packed = thinning_test._load("memberships_17ch")
    available = thinning_test._load("availability_13ch").astype(bool)
    ordering = thinning_test._load("ordering").astype(np.float32)
    shapes = thinning_test._load("registered_shapes_13ch").astype(np.float32)
    params = provenance["parameters"]
    python_result = propagate_lesion_toads_thinning(
        seg,
        packed[..., :13].astype(np.float32),
        params["object_types"],
        packed[..., 13].astype(np.float32),
        packed[..., 14].astype(np.float32),
        shapes,
        skeleton,
        available,
        ordering,
        lesion_class=params["lesion_class_index"],
        volume_ratio=params["volume_ratio"],
        limit=params["limit"],
        use_lesion_max=params["use_lesion_max"],
        source_faithful_nan=False,
    )

    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        raise RuntimeError("JDK is required for selected Java comparisons")
    with tempfile.TemporaryDirectory(prefix="p41z-iter5-to6-core-") as temp_name:
        temp = Path(temp_name)
        classes = temp / "classes"
        classes.mkdir()
        subprocess.run(
            [javac, "-d", str(classes), *(str(path) for path in (
                thinning_test.HEAP, thinning_test.LUT_JAVA,
                thinning_test.NUMERICS, thinning_test.HARNESS,
            ))],
            check=True, capture_output=True, text=True,
        )
        resource = classes / "edu/jhmi/rad/medic/structures/critical186LUT.raw.gz"
        resource.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(thinning_test.LUT_GZIP, resource)
        in_path, out_path = temp / "thinning-input.bin", temp / "thinning-java.bin"
        in_path.write_bytes(thinning_test._input_bytes(
            seg, skeleton, packed[..., :13].astype(np.float32),
            packed[..., 13].astype(np.float32), packed[..., 14].astype(np.float32),
            shapes, available, ordering,
        ))
        subprocess.run(
            [java, "-Xmx2g", "-cp", str(classes),
             "Lesion13ClassNextThinningReferenceHarness", str(in_path), str(out_path)],
            check=True, capture_output=True, text=True,
        )
        java_bytes = out_path.read_bytes()
    python_bytes = thinning_test._python_output_bytes(python_result)
    if java_bytes != python_bytes:
        raise AssertionError("iteration-6 thinning selected Java/Python serialized traces differ")

    result_paths = {
        "segmentation": THINNING_RESULT / "segmentation.nii.gz",
        "skeleton": THINNING_RESULT / "skeleton.nii.gz",
        "availability": THINNING_RESULT / "availability_13ch.nii.gz",
        "ordering": THINNING_RESULT / "ordering.nii.gz",
    }
    result_values = {
        "segmentation": python_result.run.segmentation,
        "skeleton": python_result.run.skeleton,
        "availability": python_result.run.available,
        "ordering": python_result.run.ordering,
    }
    for key, path in result_paths.items():
        if not np.array_equal(result_values[key], _load_ras(path)):
            raise AssertionError(f"fresh thinning replay differs from saved downstream {key}")

    # The next two comparators consume committed copies of the same exact
    # state. Confirm their thin output, image channels, and carried memberships
    # agree with the iteration-5 checkpoint before using those harnesses.
    case = json.loads((MEMBERSHIP_CASE / "case.json").read_text())
    files = case["fixture_files"]
    for name, metadata in files.items():
        path = MEMBERSHIP_CASE / name
        if _sha(path) != metadata["sha256"]:
            raise ValueError(f"held-fixed core fixture hash mismatch: {name}")
    np.testing.assert_array_equal(
        python_result.run.segmentation,
        _load_ras(MEMBERSHIP_CASE / case["input_files"]["segmentation"]),
    )
    np.testing.assert_array_equal(
        python_result.run.skeleton,
        _load_ras(MEMBERSHIP_CASE / "thinned_skeleton.nii.gz"),
    )
    expected_images = _load_ras(MEMBERSHIP_CASE / case["input_files"]["images"])
    expected_images_record = files[case["input_files"]["images"]]
    if expected_images_record["sha256"] != case["provenance"]["registered_recrop_t1_flair_sha256"]:
        raise ValueError("committed T1/FLAIR fixture hash does not match saved recrop provenance")
    expected_memberships = _load_ras(MEMBERSHIP_CASE / case["input_files"]["memberships"])
    np.testing.assert_array_equal(expected_memberships, packed[..., :13])
    return {
        "java_python_serialized_trace_exact": True,
        "serialized_trace_sha256": hashlib.sha256(java_bytes).hexdigest(),
        "initialization_boundary_count": len(python_result.initialization.inserted_items),
        "event_count": len(python_result.run.events),
        "downstream_state_exact": True,
        "downstream_input_crosschecks": [
            "post-thinning segmentation and skeleton",
            "registered T1/FLAIR channels",
            "iteration-5 class memberships",
        ],
        "committed_thinning_fixture_hashes": {
            key: _sha(path) for key, path in result_paths.items()
        },
        "source_selected_harness_sha256": _sha(thinning_test.HARNESS),
    }


def run() -> dict:
    thin = _run_thinning()
    stats = compare_statistics(MEMBERSHIP_CASE)
    if not stats["all_exact"]:
        raise AssertionError(f"selected statistics stage differs: {stats['difference_diagnostics']}")
    membership = compare_membership()
    if not membership["all_exact"]:
        raise AssertionError(f"selected membership stage differs: {membership['difference_diagnostics']}")
    if not membership["statistics_exact_first"]:
        raise AssertionError("membership comparator's prerequisite statistics comparison failed")
    energy_ratio = int(membership["python_energy_ratio_bits"])
    energy_ratio_value = float(np.asarray(energy_ratio, dtype=np.uint32).view(np.float32))
    report = {
        "case": "P41z synthetic saved iteration-5 to iteration-6 core replay",
        "evidence_category": "sequential selected-method Java/Python comparison; not a complete Java iteration run",
        "synthetic_only": True,
        "spatial_convention": "RAS",
        "source_checkpoint": {
            "iteration5_summary_sha256": "f03d9f0d0831022335a3296bd084e229f172ff49d86645827ec05199eb46aac0",
            "iteration5_memberships_sha256": "c6204097d018f4ac1dbcd291304878715a055a2bc2b4a8651bf4f0af30ce077a",
            "iteration5_segmentation_sha256": "f30f39161ef22a31d955db51d378e044e17acd40ebf794663e24e8570a62882d",
            "iteration5_skeleton_sha256": "2461296b1ab4d85dc114ca9df77fbf08e9ffe6131a3d9e453d7678c8999349ad",
        },
        "previous_objective": membership["previous_objective"],
        "previous_objective_provenance": membership["previous_objective_provenance"],
        "stages": {
            "iteration6_thinning": thin,
            "post_thinning_statistics": {
                "all_exact": stats["all_exact"],
                "exact_float32_bits": stats["exact_float32_bits"],
                "nan_counts": stats["nan_counts"],
                "spatial_shape_xyz": stats["spatial_shape_xyz"],
            },
            "membership_objective_and_ratio": {
                "all_exact": membership["all_exact"],
                "exact_float32_bits": membership["exact_float32_bits"],
                "nan_counts": membership["nan_counts"],
                "python_objective_bits": membership["python_objective_bits"],
                "java_objective_bits": membership["java_objective_bits"],
                "python_energy_ratio_bits": membership["python_energy_ratio_bits"],
                "java_energy_ratio_bits": membership["java_energy_ratio_bits"],
                "python_energy_ratio": energy_ratio_value,
            },
        },
        "all_selected_method_checks_exact": True,
        "scope_limits": [
            "Java references are selected-method excerpts, not a live AlgorithmLesionToads object or full composite run.",
            "The previous objective 3727.62158203125 is from the saved Python iteration-5 summary; iteration 5 itself is not independently replayed in Java.",
            "Registration, growth, and relation refresh are outside this core replay.",
        ],
    }
    RESULT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
