from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_cruise_outputs import ToadsCruiseProducts
from pymacruise.toads_dura_outputs import ToadsWrapperOutputBuffers
from pymacruise.toads_output_io import (
    save_toads_dura_wrapper_outputs,
    save_toads_outputs,
    toads_output_images,
)


SHAPE = (3, 4, 5)


def _reference(affine: np.ndarray | None = None) -> nib.Nifti1Image:
    return nib.Nifti1Image(
        np.zeros(SHAPE, dtype=np.float32),
        np.eye(4) if affine is None else affine,
    )


def _classification() -> np.ndarray:
    return np.full(SHAPE, -1, dtype=np.int8)


def _cruise_products() -> ToadsCruiseProducts:
    return ToadsCruiseProducts(
        csf=np.full(SHAPE, 0.1, dtype=np.float32),
        gray_matter=np.full(SHAPE, 0.2, dtype=np.float32),
        white_matter_fill=np.full(SHAPE, 0.3, dtype=np.float32),
        white_matter_mask=np.full(SHAPE, 1.0, dtype=np.float32),
    )


def test_hard_segmentation_outputs_unsigned_labels_and_ras_nifti(tmp_path) -> None:
    paths = save_toads_outputs(
        _reference(), _classification(), tmp_path,
        input_name="subject_t1.nii.gz", result_mode="hard_segmentation",
    )

    assert list(paths) == ["classification"]
    assert paths["classification"].name == "subject_t1_Toads_seg.nii.gz"
    image = nib.load(paths["classification"])
    assert image.get_data_dtype() == np.dtype(np.uint8)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert np.all(np.asanyarray(image.dataobj) == 255)


def test_full_segmentation_preserves_class_channels_and_optional_fields(tmp_path) -> None:
    memberships = np.stack(
        [np.full(SHAPE, value, dtype=np.float32) for value in (0.25, 0.75)],
        axis=-1,
    )
    fields = np.stack(
        [np.full(SHAPE, value, dtype=np.float32) for value in (0.9, 1.1)],
        axis=-1,
    )
    paths = save_toads_outputs(
        _reference(), _classification(), tmp_path,
        input_name="subject.nii.gz", result_mode="full_segmentation",
        memberships=memberships, inhomogeneity_field=fields,
    )

    assert paths["memberships"].name == "subject_Toads_mems.nii.gz"
    assert paths["inhomogeneity_field"].name == "subject_Toads_fields.nii.gz"
    membership_image = nib.load(paths["memberships"])
    field_image = nib.load(paths["inhomogeneity_field"])
    assert membership_image.shape == (*SHAPE, 2)
    assert field_image.shape == (*SHAPE, 2)
    assert membership_image.get_data_dtype() == np.dtype(np.float32)
    assert field_image.get_data_dtype() == np.dtype(np.float32)
    np.testing.assert_array_equal(membership_image.dataobj, memberships)
    np.testing.assert_array_equal(field_image.dataobj, fields)


@pytest.mark.parametrize("mode", ["Cruise_input", "DuraRemoval_input"])
def test_cruise_and_dura_modes_write_selected_four_products(tmp_path, mode) -> None:
    paths = save_toads_outputs(
        _reference(), _classification(), tmp_path,
        input_name="sub-01_T1", result_mode=mode,
        cruise_products=_cruise_products(),
    )

    expected_suffixes = {
        "classification": "_Toads_seg.nii.gz",
        "csf": "_csf.nii.gz",
        "gray_matter": "_gm.nii.gz",
        "white_matter_fill": "_wmfill.nii.gz",
        "white_matter_mask": "_wmmask.nii.gz",
    }
    assert set(paths) == set(expected_suffixes)
    for name, suffix in expected_suffixes.items():
        image_path = paths[name]
        assert image_path.name == f"sub-01_T1{suffix}"
        assert image_path.name.endswith(".nii.gz")
        image = nib.load(image_path)
        assert image.shape == SHAPE
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert np.all(np.asanyarray(nib.load(paths["white_matter_mask"]).dataobj) == 1.0)


def test_thalamus_mode_writes_one_float_map(tmp_path) -> None:
    thalamus = np.full((1, *SHAPE), 0.4, dtype=np.float32)
    paths = save_toads_outputs(
        _reference(), _classification(), tmp_path,
        input_name="brain.nii.gz", result_mode="thalamus_seg", thalamus=thalamus,
    )

    assert paths["thalamus"].name == "brain_thalamus.nii.gz"
    image = nib.load(paths["thalamus"])
    assert image.shape == SHAPE
    assert image.get_data_dtype() == np.dtype(np.float32)
    assert np.all(np.asanyarray(image.dataobj) == 0.4)


def test_non_ras_reference_warns_and_output_is_ras() -> None:
    affine = np.diag([-1.0, 1.0, 1.0, 1.0])
    affine[0, 3] = 2.0
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        images = toads_output_images(
            _reference(affine), _classification(), result_mode="hard_segmentation"
        )
    assert nib.aff2axcodes(images["classification"].affine) == ("R", "A", "S")


def test_rejects_wrong_mode_payloads_and_grid_shapes() -> None:
    with pytest.raises(ValueError, match="full_segmentation requires memberships"):
        toads_output_images(
            _reference(), _classification(), result_mode="full_segmentation"
        )
    with pytest.raises(ValueError, match="memberships are only exported"):
        toads_output_images(
            _reference(), _classification(), result_mode="hard_segmentation",
            memberships=np.zeros((*SHAPE, 1), dtype=np.float32),
        )
    with pytest.raises(ValueError, match="four maps matching"):
        bad = ToadsCruiseProducts(*(np.zeros((2, 4, 5), dtype=np.float32) for _ in range(4)))
        toads_output_images(
            _reference(), _classification(), result_mode="Cruise_input",
            cruise_products=bad,
        )


def test_existing_destination_is_never_overwritten(tmp_path) -> None:
    target = tmp_path / "brain_Toads_seg.nii.gz"
    target.write_bytes(b"preserve")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        save_toads_outputs(
            _reference(), _classification(), tmp_path,
            input_name="brain", result_mode="hard_segmentation",
        )
    assert target.read_bytes() == b"preserve"


def test_dura_wrapper_flat_buffers_are_reconstructed_in_mipav_x_fast_order(tmp_path) -> None:
    shape = (2, 3, 4)
    reference = nib.Nifti1Image(np.zeros(shape, dtype=np.float32), np.eye(4))
    voxel_count = int(np.prod(shape))
    labels = np.arange(voxel_count, dtype=np.int8)
    scaled_channels = tuple(
        np.arange(voxel_count, dtype=np.float32) + np.float32(channel * 100.0)
        for channel in range(4)
    )
    buffers = ToadsWrapperOutputBuffers(
        labels,
        tuple(channel / np.float32(255.0) for channel in scaled_channels),
        scaled_channels,
    )

    paths = save_toads_dura_wrapper_outputs(
        buffers, reference, tmp_path, input_name="dura-case.nii.gz"
    )

    label_image = nib.load(paths["classification"])
    assert label_image.shape == shape
    np.testing.assert_array_equal(
        label_image.dataobj, labels.reshape(shape, order="F").astype(np.uint8)
    )
    for name, buffer in zip(
        ("csf", "gray_matter", "white_matter_fill", "white_matter_mask"),
        scaled_channels,
        strict=True,
    ):
        image = nib.load(paths[name])
        np.testing.assert_array_equal(image.dataobj, buffer.reshape(shape, order="F"))
