from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_classification_update import update_lesion_classification_from_memberships
from pymacruise.lesion_post_loop_outputs import lesion_toads_post_loop_outputs
from pymacruise.lesion_state import LesionToadsFinalClassification


HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsPostLoopOutputReferenceHarness.java"
NAMES = ("SulcalCSF", "CorticalGM", "Other")


def _fixture() -> tuple[LesionToadsFinalClassification, np.ndarray, ...]:
    segmentation = np.asarray([[[0]], [[2]]], dtype=np.int8)
    memberships = np.asarray([[[[0.8, 0.1, 0.1]]], [[[0.1, 0.2, 0.7]]]], dtype=np.float32)
    lesion = np.asarray([[[0.05]], [[0.10]]], dtype=np.float32)
    black_hole = np.asarray([[[0.02]], [[0.03]]], dtype=np.float32)
    final = LesionToadsFinalClassification(
        classification=np.asarray([[[26]], [[26]]], dtype=np.int8),
        lesion_classification=np.asarray([[[1]], [[0]]], dtype=np.int16),
        changed_classification_voxels=0,
        lesion_voxels=1,
    )
    template = np.asarray([10, 20, 26], dtype=np.int8)
    return final, segmentation, memberships, lesion, black_hole, template


def _java(tmp_path: Path, mode: str, from_membership: bool, include: bool) -> dict:
    javac, java = shutil.which("javac"), shutil.which("java")
    assert javac and java, "a JDK is required for the selected-Java reference"
    classes = tmp_path / "classes"
    classes.mkdir(exist_ok=True)
    if not (classes / "LesionToadsPostLoopOutputReferenceHarness.class").exists():
        subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    result = subprocess.run(
        [java, "-cp", str(classes), "LesionToadsPostLoopOutputReferenceHarness", mode,
         str(from_membership).lower(), str(include).lower()],
        check=True, capture_output=True, text=True,
    )
    return json.loads(result.stdout)


def _run(mode: str, from_membership: bool, include: bool):
    final, segmentation, memberships, lesion, black_hole, template = _fixture()
    return lesion_toads_post_loop_outputs(
        final, memberships, lesion, black_hole, segmentation, template, NAMES, 2,
        result_mode=mode,
        output_classification_from_membership=from_membership,
        include_lesion_classification=include,
    )


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).ravel().tolist()


def _assert_products_match_java(actual, expected: dict) -> None:
    assert [product.name for product in actual.products] == [product["name"] for product in expected["products"]]
    for actual_product, java_product in zip(actual.products, expected["products"], strict=True):
        if np.issubdtype(actual_product.values.dtype, np.floating):
            assert _bits(actual_product.values) == java_product["values"]
        else:
            assert actual_product.values.ravel().tolist() == java_product["values"]


def test_dispatch_and_destination_product_order_match_java_for_all_modes_and_flags(tmp_path: Path) -> None:
    modes = (
        ("class+member", "exportMemberships"),
        ("cruise inputs", "generateCruiseInputs"),
        ("dura removal inputs", "generateDuraRemovalOutputs"),
    )
    for mode, dispatch in modes:
        for from_membership in (False, True):
            for include in (False, True):
                expected = _java(tmp_path, mode, from_membership, include)
                actual = _run(mode, from_membership, include)
                assert list(actual.source_order) == expected["order"]
                assert dispatch in actual.source_order
                assert actual.classification.ravel().tolist() == expected["classification"]
                assert actual.lesion_classification.ravel().tolist() == expected["subtype"]
                if from_membership:
                    assert actual.membership_classification.ravel().tolist() == expected["classification_mem"]
                else:
                    assert expected["classification_mem"] is None
                _assert_products_match_java(actual, expected)

                expected_names = (
                    ["memberships_plus_lesion"] if mode == "class+member" else []
                ) + ["classification"] + (
                    ["classification_from_membership"] if from_membership else []
                ) + ["lesion_subtype"] + (
                    ["cruise_inputs"] if mode == "cruise inputs" else
                    ["dura_removal_inputs"] if mode == "dura removal inputs" else []
                )
                assert [product.name for product in actual.products] == expected_names
                if mode == "cruise inputs":
                    assert actual.products[-1].values.shape == (4, 2, 1, 1)
                if mode == "dura removal inputs":
                    assert actual.products[-1].values.shape == (4, 2, 1, 1)


def test_membership_classification_snapshot_matches_existing_update_adapter() -> None:
    actual = _run("classification-only", True, False)
    final, _, memberships, lesion, _, template = _fixture()
    reference = final.classification.copy()
    update_lesion_classification_from_memberships(reference, template, memberships, lesion, 2)
    np.testing.assert_array_equal(actual.membership_classification, reference)
    assert actual.membership_classification.ravel().tolist() == [10, 26]


def test_incomplete_all_images_mode_is_rejected() -> None:
    final, segmentation, memberships, lesion, black_hole, template = _fixture()
    with pytest.raises(ValueError, match="all_images.*unsupported"):
        lesion_toads_post_loop_outputs(
            final, memberships, lesion, black_hole, segmentation, template, NAMES, 2,
            result_mode="all_images",
        )
