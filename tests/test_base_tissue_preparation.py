from __future__ import annotations

import json

import nibabel as nib
import numpy as np

from pymacruise.base_tissue_preparation import prepare_base_tissue_maps
from pymacruise.boundary import combine_multi_atlas_toads_images
from pymacruise.cli import main


def _make_inputs() -> tuple[nib.Nifti1Image, ...]:
    shape = (11, 11, 11)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    labels = np.zeros(shape, dtype=np.int16)
    labels[3:8, 3:8, 3:8] = 100  # selected MaSeg GM defines crop support
    labels[4:7, 4:7, 4:7] = 4  # selected WM, eroded by Java's 3D rule
    labels[2, 2, 2] = 31  # brain-tissue mask only, from MaSeg
    maseg = nib.Nifti1Image(labels, affine)

    # Translate source TOADS by +1 mm in RAS. Nearest-neighbor transfer to
    # MaSeg reference should move source x=0 to reference x=1.
    source_affine = affine.copy()
    source_affine[0, 3] = 1.0
    gm = np.full(shape, 0.5, dtype=np.float32)
    wm_fill = np.full(shape, 0.6, dtype=np.float32)
    toads_mask = np.full(shape, 7, dtype=np.int16)
    gm[2, 3, 3] = 0.1  # transfer maps source x=2 to reference x=3
    gm[3, 4, 4] = 0.255  # equality at threshold falls through to TOADS
    return (
        maseg,
        nib.Nifti1Image(gm, source_affine),
        nib.Nifti1Image(wm_fill, source_affine),
        nib.Nifti1Image(toads_mask, source_affine),
    )


def test_base_preparation_uses_maseg_crop_transfer_and_selected_fusion() -> None:
    maseg, toads_gm, toads_wm, toads_mask = _make_inputs()
    result = prepare_base_tissue_maps(maseg, toads_gm, toads_wm, toads_mask)

    assert result.stage_order == (
        "corticalSeparationFromMaSeg",
        "nearestNeighborTransferToMaSegReference",
        "sharedCropFromPositiveMaSegGM",
        "deriveMaSegCSFBoundary",
        "cleanBoundaryBySecondPassToadsWM",
        "combineMaSegAtlasAndTOADS",
        "cleanBoundaryByUpdatedWM",
    )
    assert result.crop_parameters.bounds == ((-2, 12), (-2, 12), (-2, 12))
    assert result.crop_parameters.cropped_shape == (15, 15, 15)
    assert result.provenance["fusion"] == {
        "operator": "CombineMAandTOADS",
        "erosion": "3D center plus six neighbors",
        "gm_confident_level": 0.001,
    }
    assert result.provenance["mask_distinction"]["fusion_mask"] == (
        "MaSeg-derived Atlas_BrainTissueMask"
    )
    assert "skull-stripping" in result.provenance["mask_distinction"]["skull_stripping_mask"]

    for image in (
        result.cropped_atlas_gm,
        result.cropped_atlas_wm,
        result.cropped_atlas_mask,
        result.cropped_atlas_brain_tissue_mask,
        result.cropped_toads_gm,
        result.cropped_toads_wm_fill,
        result.cropped_toads_mask,
        result.gray_matter,
        result.white_matter,
        result.mask,
    ):
        assert image.shape == (15, 15, 15)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(image.affine, result.cropped_atlas_gm.affine)

    # The +1 mm source-to-reference translation is nearest-neighbor sampled;
    # the source x=0 plane appears at reference x=1, then at crop x=3.
    np.testing.assert_array_equal(
        np.asarray(result.transferred_toads_gm.dataobj)[1, 3:8, 3:8],
        np.asarray(toads_gm.dataobj)[0, 3:8, 3:8],
    )
    assert np.asarray(result.cropped_atlas_brain_tissue_mask.dataobj)[4, 4, 4] == 255

    expected = combine_multi_atlas_toads_images(
        result.cropped_atlas_gm,
        result.cropped_atlas_wm,
        result.cropped_atlas_brain_tissue_mask,
        result.cropped_toads_gm,
        result.cropped_toads_wm_fill,
        result.cropped_toads_mask,
        gm_confident_level=0.001,
    )
    for actual, reference in zip(
        (result.gray_matter, result.white_matter, result.mask), expected, strict=True
    ):
        np.testing.assert_array_equal(actual.dataobj, reference.dataobj)

    # Selected CombineMAandTOADS semantics: eroded WM priority, then strict
    # low-membership atlas-GM fill, then TOADS fallback values.
    gm = np.asarray(result.gray_matter.dataobj)
    wm = np.asarray(result.white_matter.dataobj)
    mask = np.asarray(result.mask.dataobj)
    # Atlas WM is [4:7]^3, so the interior center remains after six-neighbor erosion.
    assert (gm[7, 7, 7], wm[7, 7, 7], mask[7, 7, 7]) == (0, 255, 255)
    # A selected atlas GM voxel with low positive TOADS GM is filled to 255.
    assert (gm[5, 5, 5], wm[5, 5, 5], mask[5, 5, 5]) == (255, 0, 0)
    # At exactly 0.255 the source condition is false and the TOADS value remains.
    assert gm[6, 6, 6] == np.float32(0.255)


def test_base_preparation_rejects_empty_maseg_gm_support() -> None:
    maseg, gm, wm, mask = _make_inputs()
    labels = np.zeros(maseg.shape, dtype=np.int16)
    empty = nib.Nifti1Image(labels, maseg.affine)
    try:
        prepare_base_tissue_maps(empty, gm, wm, mask)
    except ValueError as exc:
        assert "empty threshold mask" in str(exc)
    else:
        raise AssertionError("empty MaSeg GM support must not create invalid crop bounds")


def test_base_preparation_removes_selected_labels_before_transfer_and_crop() -> None:
    maseg, gm, wm, mask = _make_inputs()
    labels = np.asarray(maseg.dataobj).copy()
    labels[3, 3, 3] = 35
    maseg = nib.Nifti1Image(labels, maseg.affine)
    result = prepare_base_tissue_maps(
        maseg, gm, wm, mask, remove_cerebellum=True
    )

    assert result.stage_order == (
        "corticalSeparationFromMaSeg",
        "removeBrainstemAndCerebellumBeforeTransfer",
        "nearestNeighborTransferToMaSegReference",
        "sharedCropFromPositiveMaSegGM",
        "deriveMaSegCSFBoundary",
        "cleanBoundaryBySecondPassToadsWM",
        "combineMaSegAtlasAndTOADS",
        "cleanBoundaryByUpdatedWM",
    )
    assert result.removed_voxel_counts is not None
    assert result.removed_voxel_counts[0] >= 1
    # Label 35 is removed before TOADS transfer, and thus contributes no value
    # at the corresponding MaSeg-reference/cropped voxel.
    cropped_index = (5, 5, 5)  # original MaSeg index (3, 3, 3), crop offset -2
    assert np.asarray(result.cropped_toads_gm.dataobj)[cropped_index] == 0
    assert result.provenance["remove_cerebellum"] is True


def test_base_preparation_cli_persists_ras_maps_and_provenance(tmp_path) -> None:
    maseg, gm, wm, mask = _make_inputs()
    labels = np.asarray(maseg.dataobj).copy()
    labels[3, 3, 3] = 35
    inputs = {
        "maseg": nib.Nifti1Image(labels, maseg.affine),
        "gm": gm,
        "wm": wm,
        "mask": mask,
    }
    paths = {}
    for name, image in inputs.items():
        paths[name] = tmp_path / f"{name}.nii.gz"
        nib.save(image, paths[name])
    output_dir = tmp_path / "base-tissue"
    provenance = output_dir / "provenance.json"
    assert main([
        "base-tissue-preparation",
        "--maseg", str(paths["maseg"]),
        "--gm", str(paths["gm"]),
        "--wm-fill", str(paths["wm"]),
        "--wm-mask", str(paths["mask"]),
        "--remove-cerebellum",
        "--output-dir", str(output_dir),
        "--provenance", str(provenance),
    ]) == 0

    document = json.loads(provenance.read_text(encoding="utf-8"))
    assert document["spatial_convention"] == "RAS"
    assert document["parameters"]["remove_cerebellum"] is True
    assert document["parameters"]["removed_voxel_counts"][0] >= 1
    assert document["parameters"]["mask_distinction"]["fusion_mask"] == (
        "MaSeg-derived Atlas_BrainTissueMask"
    )
    assert set(document["outputs"]) == {
        "atlas_gm", "atlas_wm", "atlas_mask", "atlas_brain_tissue_mask",
        "toads_gm_transferred", "toads_wmfill_transferred", "toads_mask_transferred",
        "atlas_gm_cropped", "atlas_wm_cropped", "atlas_mask_cropped",
        "atlas_brain_tissue_mask_cropped", "toads_gm_cropped", "toads_wmfill_cropped",
        "toads_mask_cropped", "boundary_csf_raw", "boundary_csf_toads_cleaned",
        "boundary_csf", "fused_gm", "fused_wmfill", "fused_mask",
    }
    for name, path in document["outputs"].items():
        image = nib.load(path)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S"), name
        assert np.isfinite(np.asarray(image.dataobj)).all(), name
