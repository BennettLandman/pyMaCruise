from __future__ import annotations

import hashlib
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np

from pymacruise.lesion_post_loop_outputs import lesion_toads_post_loop_outputs
from test_lesion_nan_post_registration_branch import (
    HARNESS as BRANCH_HARNESS,
    OUTPUT_HARNESS,
    OUTPUT_HARNESS_SHA256,
    REFERENCE as BRANCH_REFERENCE,
    LABELS,
    _branch,
    _java_terminal_outputs,
    _p39ar_context,
)
from test_lesion_post_thinning_gm_membership import SPEC


ROOT = Path(__file__).parent
DISPATCH_HARNESS = ROOT / "java_reference/LesionPostRegistrationDuraDispatchReferenceHarness.java"
DISPATCH_HARNESS_SHA256 = "e130a953be23c9ecdf6cb3fce95555732050589910a0199f79854e8de23208c7"
DURA_NAMES = (
    "Background", "Sulcal-CSF", "Ventricles", "Cerebellum-GM", "Cerebrum-GM",
    "Caudate", "Thalamus", "Putamen", "Brainstem", "Cerebellum-WM",
    "OuterVentricular-WM", "InterVentricular-WM", "Frontal-WM",
)


def _write_dispatch_input(path: Path, result, java_terminal: dict, *, from_membership: bool, include: bool) -> None:
    membership = result.membership_result.membership
    shape = tuple(int(size) for size in result.segmentation.shape)
    classes = int(membership.class_memberships.shape[-1])
    header = struct.pack(
        ">7i", *shape, classes, int(SPEC["lesion_class"]), int(from_membership), int(include)
    )
    fields = (
        header,
        np.asarray(LABELS, dtype=np.int8).tobytes(),
        np.asarray(result.segmentation, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(java_terminal["classification"], dtype=np.int8).tobytes(),
        np.asarray(java_terminal["lesion_classification"], dtype=np.int16).astype(">i2").tobytes(),
        np.asarray(membership.class_memberships, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(membership.lesion_membership, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(membership.black_hole_membership, dtype=">f4").ravel(order="C").tobytes(),
    )
    path.write_bytes(b"".join(fields))


def _run_dispatch_java(classpath: Path, input_path: Path) -> dict:
    completed = subprocess.run(
        [shutil.which("java") or "java", "-cp", str(classpath),
         "LesionPostRegistrationDuraDispatchReferenceHarness", str(input_path)],
        check=True, capture_output=True, text=True,
    )
    import json
    return json.loads(completed.stdout)


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).ravel().tolist()


def test_p39ar_terminal_state_dispatches_dura_outputs_in_java_destination_order(tmp_path: Path) -> None:
    assert hashlib.sha256(BRANCH_HARNESS.read_bytes()).hexdigest() == BRANCH_REFERENCE["harness_sha256"]
    assert hashlib.sha256(OUTPUT_HARNESS.read_bytes()).hexdigest() == OUTPUT_HARNESS_SHA256
    assert hashlib.sha256(DISPATCH_HARNESS.read_bytes()).hexdigest() == DISPATCH_HARNESS_SHA256
    javac = shutil.which("javac")
    assert javac is not None and shutil.which("java") is not None, "a JDK is required for selected-Java references"
    subprocess.run([javac, "-d", str(tmp_path), str(OUTPUT_HARNESS)], check=True, capture_output=True, text=True)
    subprocess.run([javac, "-d", str(tmp_path), str(DISPATCH_HARNESS)], check=True, capture_output=True, text=True)

    frames, update, registration, thinning = _p39ar_context()
    assert frames.shape[:3] == (17, 17, 17)
    assert len(DURA_NAMES) == len(SPEC["object_types"]) == 13
    assert len(BRANCH_HARNESS.read_bytes()) > 0
    result = _branch(
        update, registration, thinning,
        previous=float(registration.energy_ratio.previous_energy), iteration=3,
        n_iter_max=-3, max_distance=0.1, source_ratio=registration.energy_ratio,
    )
    assert result.stop and not result.run_growing and result.final_classification is not None
    assert result.stage_order[-2:] == ("updateClassificationFromLabels", "updateLesionClassifcation")

    terminal_java = _java_terminal_outputs(
        result, result.membership_result, classpath=tmp_path,
        output_path=tmp_path / "p39ar-terminal-state.bin",
    )
    np.testing.assert_array_equal(
        result.final_classification.classification.ravel(),
        np.asarray(terminal_java["classification"], dtype=np.int16),
    )
    np.testing.assert_array_equal(
        result.final_classification.lesion_classification.ravel(),
        np.asarray(terminal_java["lesion_classification"], dtype=np.int16),
    )

    for from_membership in (False, True):
        for include in (False, True):
            input_path = tmp_path / f"dispatch-{int(from_membership)}-{int(include)}.bin"
            _write_dispatch_input(input_path, result, terminal_java,
                                  from_membership=from_membership, include=include)
            expected = _run_dispatch_java(tmp_path, input_path)
            actual = lesion_toads_post_loop_outputs(
                result.final_classification,
                result.membership_result.membership.class_memberships,
                result.membership_result.membership.lesion_membership,
                result.membership_result.membership.black_hole_membership,
                result.segmentation,
                LABELS,
                DURA_NAMES,
                int(SPEC["lesion_class"]),
                result_mode="dura removal inputs",
                output_classification_from_membership=from_membership,
                include_lesion_classification=include,
            )
            expected_products = expected["products"]
            expected_names = [product["name"] for product in expected_products]
            adapter_names = [product.name for product in actual.products]
            assert adapter_names == (
                ["classification"]
                + (["classification_from_membership"] if from_membership else [])
                + ["lesion_subtype", "dura_removal_inputs"]
            )
            expanded_names = adapter_names[:-1] + [
                "csf", "gray_matter", "white_matter_fill", "white_matter_mask"
            ]
            assert expanded_names == expected_names
            assert list(actual.source_order) == [
                "cleanUp", "exportClassification", "exportLesionClassification",
                *( ["updateClassificationFromMemberships", "getClassification"] if from_membership else [] ),
                "generateDuraRemovalOutputs",
                *( ["includeLesionnClassification"] if include else [] ),
            ]
            for product, java_product in zip(actual.products[:-1], expected_products[:-4], strict=True):
                assert product.values.ravel().tolist() == java_product["values"]
            assert actual.products[0].values.ravel().tolist() == expected_products[0]["values"]
            if from_membership:
                assert actual.products[1].values.ravel().tolist() == expected_products[1]["values"]
            assert actual.products[-2].values.ravel().tolist() == expected_products[-5]["values"]
            expected_channel_bits = [value for channel in expected_products[-4:] for value in channel["values"]]
            assert _bits(actual.products[-1].values) == expected_channel_bits
