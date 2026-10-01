from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_output_io import (
    lesion_toads_output_images,
    save_lesion_toads_outputs,
)
from pymacruise.lesion_post_loop_outputs import (
    LesionToadsOutputProduct,
    LesionToadsPostLoopOutputs,
)


def _outputs(mode: str = "dura removal inputs") -> LesionToadsPostLoopOutputs:
    shape = (3, 4, 5)
    channels = np.stack(
        [np.full(shape, value, dtype=np.float32) for value in (0.1, 0.2, 0.3, 0.4)]
    )
    products = (
        LesionToadsOutputProduct("classification", np.full(shape, -1, dtype=np.int8)),
        LesionToadsOutputProduct("lesion_subtype", np.full(shape, 32000, dtype=np.int16)),
        LesionToadsOutputProduct("dura_removal_inputs", channels),
    )
    return LesionToadsPostLoopOutputs(
        classification=products[0].values,
        membership_classification=None,
        lesion_classification=products[1].values,
        result_mode=mode,
        products=products,
        source_order=("classification", "lesion_subtype", "dura_removal_inputs"),
    )


def test_builds_typed_output_images_on_ras_reference_grid() -> None:
    shape = (3, 4, 5)
    reference = nib.Nifti1Image(np.zeros(shape, dtype=np.float32), np.diag([2.0, 3.0, 4.0, 1.0]))

    images = lesion_toads_output_images(_outputs(), reference)

    assert set(images) == {"classification", "lesion_subtype", "csf", "gm", "wmFill", "wmMask"}
    for image in images.values():
        assert image.shape == shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(image.affine, reference.affine)
    assert images["classification"].get_data_dtype() == np.dtype(np.uint8)
    assert images["lesion_subtype"].get_data_dtype() == np.dtype(np.int16)
    assert images["gm"].get_data_dtype() == np.dtype(np.float32)
    assert np.all(np.asanyarray(images["classification"].dataobj) == 255)
    assert np.all(np.asanyarray(images["lesion_subtype"].dataobj) == 32000)
    assert np.all(np.asanyarray(images["wmMask"].dataobj) == np.float32(0.4))


def test_class_and_member_mode_preserves_fourth_channel_dimension() -> None:
    shape = (3, 4, 5)
    stack = np.arange(np.prod(shape) * 3, dtype=np.float32).reshape((*shape, 3))
    outputs = LesionToadsPostLoopOutputs(
        classification=np.zeros(shape, dtype=np.int8),
        membership_classification=np.ones(shape, dtype=np.int8),
        lesion_classification=np.zeros(shape, dtype=np.int16),
        result_mode="class+member",
        products=(
            LesionToadsOutputProduct("memberships_plus_lesion", stack),
            LesionToadsOutputProduct("classification", np.zeros(shape, dtype=np.int8)),
            LesionToadsOutputProduct("classification_from_membership", np.ones(shape, dtype=np.int8)),
            LesionToadsOutputProduct("lesion_subtype", np.zeros(shape, dtype=np.int16)),
        ),
        source_order=("exportMemberships", "classification"),
    )
    reference = nib.Nifti1Image(np.zeros(shape, dtype=np.float32), np.eye(4))

    image = lesion_toads_output_images(outputs, reference)["memberships_plus_lesion"]

    assert image.shape == (*shape, 3)
    assert image.get_data_dtype() == np.dtype(np.float32)
    np.testing.assert_array_equal(image.dataobj, stack)


@pytest.mark.parametrize(
    ("field", "suffix", "expected_shape"),
    [
        (np.full((3, 4, 5), 1.0, dtype=np.float64), "_lToads_field.nii.gz", (3, 4, 5)),
        (np.stack([np.full((3, 4, 5), value) for value in (0.8, 1.2)], axis=-1), "_lToads_fields.nii.gz", (3, 4, 5, 2)),
    ],
)
def test_exports_optional_single_or_multichannel_restored_field(
    tmp_path, field, suffix, expected_shape
) -> None:
    reference = nib.Nifti1Image(np.zeros((3, 4, 5), dtype=np.float32), np.eye(4))

    images = lesion_toads_output_images(
        _outputs(), reference, inhomogeneity_field=field
    )
    paths = save_lesion_toads_outputs(
        _outputs(), reference, tmp_path, input_name="case.nii.gz",
        inhomogeneity_field=field,
    )

    image = images["inhomogeneity_field"]
    assert image.shape == expected_shape
    assert image.get_data_dtype() == np.dtype(np.float32)
    assert paths["inhomogeneity_field"].name == f"case{suffix}"
    restored = nib.load(paths["inhomogeneity_field"])
    assert restored.shape == expected_shape
    assert restored.get_data_dtype() == np.dtype(np.float32)
    np.testing.assert_array_equal(restored.dataobj, np.asarray(field, dtype=np.float32))


def test_saves_only_nii_gz_with_ras_geometry_and_selected_suffixes(tmp_path) -> None:
    shape = (3, 4, 5)
    reference = nib.Nifti1Image(np.zeros(shape, dtype=np.float32), np.diag([1.2, 1.4, 2.0, 1.0]))

    paths = save_lesion_toads_outputs(
        _outputs(), reference, tmp_path, input_name="case-01_t1.nii.gz"
    )

    assert paths["classification"].name == "case-01_t1_lToads_class.nii.gz"
    assert paths["lesion_subtype"].name == "case-01_t1_lToads_lesions_hard.nii.gz"
    assert paths["wmFill"].name == "case-01_t1_lToads_wmFill.nii.gz"
    for name, path in paths.items():
        assert path.exists()
        assert path.name.endswith(".nii.gz")
        restored = nib.load(path)
        assert nib.aff2axcodes(restored.affine) == ("R", "A", "S")
        np.testing.assert_allclose(restored.affine, reference.affine, rtol=0, atol=1e-6)
        assert restored.shape == shape
        expected_dtype = np.uint8 if name == "classification" else (
            np.int16 if name == "lesion_subtype" else np.float32
        )
        assert restored.get_data_dtype() == np.dtype(expected_dtype)


def test_warns_and_reorients_non_ras_reference_before_output_creation() -> None:
    shape = (3, 4, 5)
    affine = np.diag([-2.0, 3.0, 4.0, 1.0])
    affine[0, 3] = 4.0
    reference = nib.Nifti1Image(np.zeros(shape, dtype=np.float32), affine)

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        images = lesion_toads_output_images(_outputs(), reference)

    assert nib.aff2axcodes(images["classification"].affine) == ("R", "A", "S")
    assert images["classification"].shape == shape
    np.testing.assert_array_equal(images["classification"].affine, np.diag([2.0, 3.0, 4.0, 1.0]))


def test_rejects_existing_output_instead_of_overwriting(tmp_path) -> None:
    reference = nib.Nifti1Image(np.zeros((3, 4, 5), dtype=np.float32), np.eye(4))
    target = tmp_path / "case_lToads_class.nii.gz"
    target.write_bytes(b"preserve")

    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        save_lesion_toads_outputs(_outputs(), reference, tmp_path, input_name="case.nii.gz")
    assert target.read_bytes() == b"preserve"


def test_rejects_output_arrays_not_matching_the_t1_grid() -> None:
    outputs = _outputs()
    bad_product = LesionToadsOutputProduct("classification", np.zeros((2, 4, 5), dtype=np.int8))
    invalid = LesionToadsPostLoopOutputs(
        classification=bad_product.values,
        membership_classification=None,
        lesion_classification=outputs.lesion_classification,
        result_mode=outputs.result_mode,
        products=(bad_product, *outputs.products[1:]),
        source_order=outputs.source_order,
    )
    reference = nib.Nifti1Image(np.zeros((3, 4, 5), dtype=np.float32), np.eye(4))

    with pytest.raises(ValueError, match="classification must match"):
        lesion_toads_output_images(invalid, reference)


def test_rejects_empty_or_wrong_grid_optional_field() -> None:
    outputs = _outputs()
    reference = nib.Nifti1Image(np.zeros((3, 4, 5), dtype=np.float32), np.eye(4))

    for field in (np.zeros((2, 4, 5)), np.zeros((3, 4, 5, 0))):
        with pytest.raises(ValueError, match="inhomogeneity_field must be"):
            lesion_toads_output_images(
                outputs, reference, inhomogeneity_field=field
            )
