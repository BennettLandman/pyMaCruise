from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib


ROOT = Path(__file__).parents[1]
ATLAS_ROOT = ROOT / "data" / "atlases"


def test_bundled_atlases_are_hashed_ras_nifti_and_source_matched() -> None:
    audit = json.loads((ATLAS_ROOT / "CONVERSION_AUDIT.json").read_text())
    checksum_entries = {
        line.split("  ", 1)[1]: line.split("  ", 1)[0]
        for line in (ATLAS_ROOT / "SHA256SUMS").read_text().splitlines()
    }
    assert len(audit["images"]) == 35
    assert len(list(ATLAS_ROOT.rglob("*.nii.gz"))) == 35
    assert not list(ATLAS_ROOT.rglob("*.raw"))
    assert not list(ATLAS_ROOT.rglob("*.img"))
    assert not list(ATLAS_ROOT.rglob("*.hdr"))

    for row in audit["images"]:
        image_path = ATLAS_ROOT / row["nifti_output"]
        relative_path = image_path.relative_to(ATLAS_ROOT).as_posix()
        image = nib.load(str(image_path))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert row["output_axis_codes"] == ["R", "A", "S"]
        assert row["source_voxel_values_exact_after_ras_reorientation"] is True
        assert checksum_entries[relative_path] == hashlib.sha256(image_path.read_bytes()).hexdigest()


def test_atlas_descriptors_reference_only_compressed_nifti() -> None:
    for descriptor in ATLAS_ROOT.rglob("*.txt"):
        content = descriptor.read_text()
        assert ".raw" not in content
        assert ".img" not in content
        assert not content or ".nii.gz" in content or descriptor.name.startswith("oasis-")
