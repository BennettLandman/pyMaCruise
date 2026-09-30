from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_classification import generate_toads_classification


DATA = Path(__file__).parent / "data/synthetic/toads_generate_classification"
REFERENCE = DATA / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference/ToadsGenerateClassificationReferenceHarness.java"
SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
HARNESS_SHA = "e6121d3828eeda1be15f35440d2ea1c8b8ccda8d1284cc84f84ba37db3a97f3c"
FIXTURE_SHA = "bd5005b403e2cef6f4152c213690721736b8a799aeda5398834bd19eb5153d60"


def _load_inputs() -> tuple[np.ndarray, np.ndarray]:
    membership_image = nib.load(str(DATA / "memberships.nii.gz"))
    best_image = nib.load(str(DATA / "best_indices.nii.gz"))
    assert nib.aff2axcodes(membership_image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(best_image.affine) == ("R", "A", "S")
    return np.asarray(membership_image.dataobj, dtype=np.float32), np.asarray(best_image.dataobj, dtype=np.int16)


def test_java_excerpt_and_reference_fixture_are_pinned() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == FIXTURE_SHA
    assert json.loads(REFERENCE.read_text())["source_sha256"] == SOURCE_SHA


def test_exact_and_approximate_output_match_selected_java_excerpt() -> None:
    reference = json.loads(REFERENCE.read_text())
    memberships, best = _load_inputs()
    template = np.asarray(reference["template_labels"], dtype=np.int8)
    fixture_memberships = np.transpose(memberships, (2, 1, 0, 3)).reshape(-1)
    encoded_memberships = [
        "NaN" if np.isnan(value) else
        "Infinity" if np.isposinf(value) else
        "-Infinity" if np.isneginf(value) else
        float(format(float(value), ".7g"))
        for value in fixture_memberships
    ]
    assert encoded_memberships == reference["memberships_x_fastest_class_inner"]
    assert best.ravel(order="F").tolist() == reference["best_indices_x_fastest"]
    exact = generate_toads_classification(template, memberships)
    approx = generate_toads_classification(
        template, memberships, use_approx=True, best_class_indices=best
    )
    assert exact.dtype == np.int8 and approx.dtype == np.int8
    assert -128 in approx and -1 in approx
    np.testing.assert_array_equal(exact.ravel(order="F"), reference["exact_x_fastest"])
    np.testing.assert_array_equal(approx.ravel(order="F"), reference["approx_x_fastest"])
    assert exact[0, 0, 0] == 0  # all scores are zero
    assert exact[1, 0, 0] == 0  # all scores are negative
    assert exact[2, 0, 0] == 0  # NaN-only scores never compare greater
    assert exact[0, 1, 0] == 42  # positive infinity wins
    assert exact[1, 1, 0] == -1  # a finite class wins the first tie
    assert exact[2, 1, 0] == 127
    assert exact[1, 0, 1] == -1  # first positive infinity wins
    assert exact[2, 0, 1] == -1  # ascending class scan keeps class 1 on a tie
    assert approx[0, 0, 0] == -1  # approximate path maps best index to signed label
