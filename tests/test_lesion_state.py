from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_state import initialize_lesion_toads_segmentation_state


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "lesion_segmentation_state"
CLASS_LABELS = (1, 2, 3)
OBJECT_TYPES = ("obj", "lesion", "optimized_blackHole")


def _fixture_inputs() -> tuple[np.ndarray, np.ndarray]:
    t1 = nib.load(FIXTURE / "t1.nii.gz")
    flair = nib.load(FIXTURE / "flair.nii.gz")
    labels_image = nib.load(FIXTURE / "template_labels.nii.gz")
    for image in (t1, flair, labels_image):
        np.testing.assert_allclose(image.affine, np.eye(4))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    channels = np.stack((np.asanyarray(t1.dataobj), np.asanyarray(flair.dataobj)))
    labels = np.asanyarray(labels_image.dataobj)
    return channels.astype(np.float32), labels


def _flat(array: np.ndarray) -> list[float | int | bool]:
    return array.ravel(order="C").tolist()


def test_constructor_state_matches_selected_java_formula_fixture() -> None:
    channels, labels = _fixture_inputs()
    reference = json.loads((FIXTURE / "java_reference.json").read_text())
    np.testing.assert_array_equal(
        labels.ravel(order="C"), np.asarray(reference["labels"], dtype=np.int16)
    )
    np.testing.assert_array_equal(
        np.stack([channel.ravel(order="C") for channel in channels]).ravel(),
        np.asarray(reference["images"], dtype=np.float32),
    )

    state = initialize_lesion_toads_segmentation_state(
        channels, labels, CLASS_LABELS, OBJECT_TYPES
    )

    for name in ("segmentation", "skeleton", "relation_map"):
        expected_dtype = np.int16 if name == "relation_map" else np.int8
        np.testing.assert_array_equal(
            _flat(getattr(state, name)), np.asarray(reference[name], dtype=expected_dtype)
        )
    np.testing.assert_array_equal(
        _flat(state.memberships), np.asarray(reference["memberships"], dtype=np.float32)
    )
    np.testing.assert_array_equal(_flat(state.available), reference["available"])
    np.testing.assert_array_equal(state.class_volume, reference["class_volume"])
    expected_centroid = np.asarray(reference["image_centroid"], dtype=object).reshape(2, 3)
    for channel, class_index in ((0, 0), (0, 1), (1, 0), (1, 1)):
        assert state.image_centroid[channel, class_index] == np.float32(
            expected_centroid[channel, class_index]
        )
    assert np.isnan(state.image_centroid[:, 2]).all()
    assert reference["absent_class_centroid_is_nan"] is True
    assert state.lesion_class == reference["lesion_class"] == 1
    np.testing.assert_array_equal(state.lesion_similarity, reference["lesion_similarity"])
    np.testing.assert_array_equal(
        state.black_hole_similarity, reference["black_hole_similarity"]
    )
    np.testing.assert_array_equal(
        state.class_similarity.ravel(), reference["class_similarity"]
    )
    assert np.all(state.intensity_stddev == np.float32(1.0))
    assert np.all(state.ordering == np.float32(0.0))
    assert np.all(state.lesion_classification == 0)
    assert np.all(state.lesion_classification_membership == 0)
    assert np.all(state.lesion_membership == np.float32(0.0))
    assert np.all(state.lesion_relation == np.float32(1.0))
    assert np.all(state.outlier_membership == np.float32(0.5))
    assert np.all(state.black_hole_membership == np.float32(0.0))
    assert np.all(state.ventricle_distance == reference["ventricle_sentinel"])
    assert np.all(state.gray_matter_distance == reference["gm_sentinel"])
    assert np.all(state.brainstem_distance == reference["brainstem_sentinel"])
    assert np.all(state.csf_distance == 0)
    for mask in (
        state.ventricle_mask,
        state.gray_matter_mask,
        state.brainstem_mask,
        state.csf_mask,
    ):
        assert not np.any(mask)
    np.testing.assert_array_equal(state.covariance, np.broadcast_to(np.eye(2), (3, 2, 2)))
    np.testing.assert_array_equal(state.lesion_covariance, np.eye(2))
    assert np.all(state.prior_centroid == 0)
    assert np.all(state.lesion_centroid == 0)
    assert np.all(state.black_hole_centroid == 0)
    assert np.all(state.lesion_stddev == 1)
    assert np.all(state.black_hole_stddev == 1)
    assert np.all(state.image_high == 1)
    assert np.all(state.image_scale == 1)
    assert np.all(state.channel_variance == 1)
    np.testing.assert_array_equal(state.channel_weight, np.full(2, 0.5, np.float32))
    np.testing.assert_array_equal(state.lesion_weight, np.full(2, 0.5, np.float32))


def test_no_lesion_classes_keep_default_similarity_and_class_sentinel() -> None:
    channels, labels = _fixture_inputs()
    reference = json.loads((FIXTURE / "java_reference.json").read_text())
    state = initialize_lesion_toads_segmentation_state(
        channels, labels, CLASS_LABELS, ("obj", "mask", "out")
    )
    assert state.lesion_class == reference["no_lesion_class"] == -1
    np.testing.assert_array_equal(
        state.lesion_similarity, reference["no_lesion_similarity"]
    )
    np.testing.assert_array_equal(
        state.black_hole_similarity, reference["no_lesion_black_hole_similarity"]
    )


@pytest.mark.parametrize(
    ("images", "labels", "class_labels", "types", "message"),
    [
        (np.zeros((1, 2, 2, 2)), np.ones((2, 2, 2), dtype=np.int16), (0,), ("obj",), "absent"),
        (np.zeros((1, 2, 2, 1)), np.zeros((2, 2, 2), dtype=np.int16), (0,), ("obj",), "grids must match"),
        (np.zeros((1, 2, 2, 2)), np.zeros((2, 2, 2), dtype=np.float32), (0,), ("obj",), "must be integers"),
        (np.zeros((1, 2, 2, 2)), np.zeros((2, 2, 2), dtype=np.int16), (0,), ("obj", "mask"), "one object type"),
    ],
)
def test_constructor_state_rejects_invalid_contract(
    images: np.ndarray,
    labels: np.ndarray,
    class_labels: tuple[int, ...],
    types: tuple[str, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        initialize_lesion_toads_segmentation_state(images, labels, class_labels, types)
