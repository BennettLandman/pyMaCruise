from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_state import lesion_toads_final_classification


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_final_classification"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / f"{name}.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj)


def test_final_lesion_classification_matches_selected_java_methods() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert REFERENCE["shape"] == [3, 2, 2]
    segmentation = _load("segmentation")
    previous = _load("previous_classification")
    memberships = _load("memberships")
    lesion = _load("lesion_membership")
    black_hole = _load("black_hole_membership")
    source_inputs = REFERENCE["inputs"]
    np.testing.assert_array_equal(segmentation.ravel(), source_inputs["segmentation"])
    np.testing.assert_array_equal(
        previous.ravel(), source_inputs["previous_classification"]
    )
    np.testing.assert_array_equal(
        memberships.ravel(), np.asarray(source_inputs["memberships"], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        lesion.ravel(), np.asarray(source_inputs["lesion_membership"], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        black_hole.ravel(),
        np.asarray(source_inputs["black_hole_membership"], dtype=np.float32),
    )

    for reference in REFERENCE["cases"]:
        result = lesion_toads_final_classification(
            segmentation,
            REFERENCE["template_labels"],
            previous,
            memberships,
            lesion,
            black_hole,
            lesion_class=1,
            has_lesion=reference["has_lesion"],
        )
        np.testing.assert_array_equal(
            result.classification.ravel(), reference["classification"]
        )
        np.testing.assert_array_equal(
            result.lesion_classification.ravel(), reference["lesion_classification"]
        )
        assert result.changed_classification_voxels == reference["class_diff"]
        assert result.lesion_voxels == reference["lesion_volume"]


@pytest.mark.parametrize(
    "changes, message",
    [
        ({"lesion_class": 3}, "outside the membership class axis"),
        ({"has_lesion": 1}, "has_lesion must be boolean"),
    ],
)
def test_final_lesion_classification_rejects_invalid_arguments(changes, message) -> None:
    kwargs = {"lesion_class": 1, "has_lesion": True}
    kwargs.update(changes)
    with pytest.raises(ValueError, match=message):
        lesion_toads_final_classification(
            _load("segmentation"),
            REFERENCE["template_labels"],
            _load("previous_classification"),
            _load("memberships"),
            _load("lesion_membership"),
            _load("black_hole_membership"),
            **kwargs,
        )


def test_final_lesion_classification_rejects_mismatched_auxiliary_grid() -> None:
    with pytest.raises(ValueError, match="must match segmentation"):
        lesion_toads_final_classification(
            _load("segmentation"),
            REFERENCE["template_labels"],
            _load("previous_classification"),
            _load("memberships"),
            _load("lesion_membership")[..., :1],
            _load("black_hole_membership"),
            lesion_class=1,
        )
