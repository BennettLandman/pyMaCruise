from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np
import pytest

from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_images


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "tests" / "data" / "synthetic" / "tgdm"
PRODUCTS = (
    "initial_level_set",
    "tissue_sign",
    "initial_inner",
    "inner",
    "central",
    "outer",
)


def _save_lps_equivalent(image: nib.spatialimages.SpatialImage, path: Path) -> None:
    """Store the same physical samples in flipped-x, non-RAS voxel order."""
    data = np.asarray(image.dataobj)
    flip = np.eye(4)
    flip[0, 0] = -1.0
    flip[0, 3] = data.shape[0] - 1
    nib.save(nib.Nifti1Image(data[::-1, ...].copy(), image.affine @ flip), path)


def test_tgdm_cli_matches_api_for_all_products_optional_inputs_and_ras_provenance(
    tmp_path: Path,
) -> None:
    images = {
        name: nib.load(FIXTURE / f"{name}.nii.gz")
        for name in (
            "initial_membership",
            "white_matter",
            "gray_matter",
            "ace_gray_matter",
            "ace_skeleton",
        )
    }
    reference = images["initial_membership"]
    grid = np.indices(reference.shape, dtype=np.float32)
    center = np.asarray(reference.shape, dtype=np.float32).reshape((3, 1, 1, 1)) / 2.0
    radius = np.sqrt(np.sum((grid - center) ** 2, axis=0))
    brain_mask = nib.Nifti1Image((radius <= 10).astype(np.uint8), reference.affine)
    gvf = nib.Nifti1Image(
        np.zeros((*reference.shape, 3), dtype=np.float32), reference.affine
    )

    # The optional ACE skeleton is stored with a non-RAS voxel axis order.
    # Canonicalizing it should recover the same physical image while warning.
    lps_skeleton_path = tmp_path / "ace-skeleton-lps.nii.gz"
    _save_lps_equivalent(images["ace_skeleton"], lps_skeleton_path)
    brain_mask_path = tmp_path / "brain-mask.nii.gz"
    gvf_path = tmp_path / "gvf.nii.gz"
    nib.save(brain_mask, brain_mask_path)
    nib.save(gvf, gvf_path)

    parameters = TGDMBaselineParameters()
    with pytest.warns(UserWarning, match="reorienting the image to RAS") as caught:
        expected = tgdm_nested_images(
            images["initial_membership"],
            images["white_matter"],
            images["gray_matter"],
            images["ace_gray_matter"],
            brain_mask=brain_mask,
            ace_skeleton=nib.load(lps_skeleton_path),
            gvf=gvf,
            parameters=parameters,
        )
    assert len(caught) == 1

    output_dir = tmp_path / "tgdm-output"
    provenance_path = tmp_path / "tgdm.json"
    command = [
        sys.executable,
        "-m",
        "pymacruise.cli",
        "tgdm",
        "--initial",
        str(FIXTURE / "initial_membership.nii.gz"),
        "--white",
        str(FIXTURE / "white_matter.nii.gz"),
        "--gray",
        str(FIXTURE / "gray_matter.nii.gz"),
        "--ace-gray",
        str(FIXTURE / "ace_gray_matter.nii.gz"),
        "--brain-mask",
        str(brain_mask_path),
        "--ace-skeleton",
        str(lps_skeleton_path),
        "--gvf",
        str(gvf_path),
        "--output-dir",
        str(output_dir),
        "--provenance",
        str(provenance_path),
    ]
    completed = subprocess.run(
        command, cwd=REPO, capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr
    assert "ACE skeleton input orientation" in completed.stderr
    assert "reorienting the image to RAS" in completed.stderr

    for name in PRODUCTS:
        actual = nib.load(output_dir / f"{name}.nii.gz")
        expected_image = getattr(expected, name)
        assert nib.aff2axcodes(actual.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(actual.affine, expected_image.affine)
        np.testing.assert_array_equal(actual.dataobj, expected_image.dataobj)

    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["schema_version"] == 1
    assert provenance["command"] == "tgdm"
    assert provenance["software"]["name"] == "pymacruise"
    assert provenance["software"]["version"]
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["inputs"] == {
        "initial_membership": str((FIXTURE / "initial_membership.nii.gz").resolve()),
        "white_matter": str((FIXTURE / "white_matter.nii.gz").resolve()),
        "gray_matter": str((FIXTURE / "gray_matter.nii.gz").resolve()),
        "ace_gray_matter": str((FIXTURE / "ace_gray_matter.nii.gz").resolve()),
        "brain_mask": str(brain_mask_path.resolve()),
        "ace_skeleton": str(lps_skeleton_path.resolve()),
        "gvf": str(gvf_path.resolve()),
    }
    assert provenance["parameters"] == asdict(parameters)
    assert provenance["outputs"] == {
        name: str((output_dir / f"{name}.nii.gz").resolve()) for name in PRODUCTS
    }
