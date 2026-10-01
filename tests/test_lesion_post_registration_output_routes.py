from __future__ import annotations

import hashlib
import json
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np

from pymacruise.lesion_post_loop_outputs import lesion_toads_post_loop_outputs
from test_lesion_nan_post_registration_branch import (
    OUTPUT_HARNESS,
    OUTPUT_HARNESS_SHA256,
    _branch,
    _java_terminal_outputs,
    _p39ar_context,
)
from test_lesion_post_thinning_gm_membership import SPEC


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference/LesionPostRegistrationOutputRoutesReferenceHarness.java"
HARNESS_SHA256 = "ca1dd209040a84c7f0091fdf3f2f33e6ff3d9d4e3ab3ce6a97046035abf0cab2"
NAMES = (
    "Background", "Sulcal-CSF", "Ventricles", "Cerebellum-GM", "Cerebrum-GM",
    "Caudate", "Thalamus", "Putamen", "Brainstem", "Cerebellum-WM",
    "OuterVentricular-WM", "InterVentricular-WM", "Frontal-WM",
)
LABELS = np.arange(10, 10 + len(SPEC["object_types"]), dtype=np.int8)


def _write_input(path: Path, result, terminal: dict, *, from_membership: bool, include: bool) -> None:
    membership = result.membership_result.membership
    shape = tuple(int(size) for size in result.segmentation.shape)
    classes = int(membership.class_memberships.shape[-1])
    parts = (
        struct.pack(">7i", *shape, classes, int(SPEC["lesion_class"]), int(from_membership), int(include)),
        LABELS.tobytes(),
        np.asarray(result.segmentation, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(terminal["classification"], dtype=np.int8).tobytes(),
        np.asarray(terminal["lesion_classification"], dtype=np.int16).astype(">i2").tobytes(),
        np.asarray(membership.class_memberships, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(membership.lesion_membership, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(membership.black_hole_membership, dtype=">f4").ravel(order="C").tobytes(),
    )
    path.write_bytes(b"".join(parts))


def _java(classpath: Path, path: Path, mode: str) -> dict:
    run = subprocess.run(
        [shutil.which("java") or "java", "-cp", str(classpath),
         "LesionPostRegistrationOutputRoutesReferenceHarness", str(path), mode],
        check=True, capture_output=True, text=True,
    )
    return json.loads(run.stdout)


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).ravel().tolist()


def test_p39ar_terminal_state_composes_class_member_and_cruise_routes(tmp_path: Path) -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(OUTPUT_HARNESS.read_bytes()).hexdigest() == OUTPUT_HARNESS_SHA256
    javac = shutil.which("javac")
    assert javac is not None and shutil.which("java") is not None, "a JDK is required for selected-Java references"
    subprocess.run(
        [javac, "-d", str(tmp_path), str(OUTPUT_HARNESS), str(HARNESS)],
        check=True, capture_output=True, text=True,
    )

    frames, update, registration, thinning = _p39ar_context()
    assert frames.shape[:3] == (17, 17, 17)
    assert len(NAMES) == len(SPEC["object_types"]) == 13
    branch = _branch(
        update, registration, thinning,
        previous=float(registration.energy_ratio.previous_energy), iteration=3,
        n_iter_max=-3, max_distance=0.1, source_ratio=registration.energy_ratio,
    )
    assert branch.stop and not branch.run_growing and branch.final_classification is not None
    terminal = _java_terminal_outputs(
        branch, branch.membership_result, classpath=tmp_path,
        output_path=tmp_path / "terminal.bin",
    )

    for mode in ("class+member", "cruise inputs"):
        for from_membership in (False, True):
            for include in (False, True):
                case_path = tmp_path / f"{mode.replace(' ', '-')}-{int(from_membership)}-{int(include)}.bin"
                _write_input(case_path, branch, terminal,
                             from_membership=from_membership, include=include)
                expected = _java(tmp_path, case_path, mode)
                actual = lesion_toads_post_loop_outputs(
                    branch.final_classification,
                    branch.membership_result.membership.class_memberships,
                    branch.membership_result.membership.lesion_membership,
                    branch.membership_result.membership.black_hole_membership,
                    branch.segmentation,
                    LABELS,
                    NAMES,
                    int(SPEC["lesion_class"]),
                    result_mode=mode,
                    output_classification_from_membership=from_membership,
                    include_lesion_classification=include,
                )
                assert list(actual.source_order) == expected["order"]
                actual_products = actual.products
                java_products = expected["products"]
                assert [p.name for p in actual_products] == [p["name"] for p in java_products]
                for product, java_product in zip(actual_products, java_products, strict=True):
                    if np.issubdtype(product.values.dtype, np.floating):
                        assert _bits(product.values) == java_product["values"]
                    else:
                        assert product.values.ravel().tolist() == java_product["values"]
