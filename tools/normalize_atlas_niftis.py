"""Normalize user-converted atlas NIfTIs to RAS and their source folders.

This one-time repository maintenance utility keeps encoded voxel values and
physical locations while permuting/flipping voxel axes only. It records both
the user-provided NIfTI affine and the final RAS affine in a JSON inventory.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import (
    axcodes2ornt,
    apply_orientation,
    inv_ornt_aff,
    io_orientation,
    ornt_transform,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _affine_list(affine: np.ndarray) -> list[list[float]]:
    return np.asarray(affine, dtype=np.float64).tolist()


def normalize_atlases(atlas_root: Path) -> list[dict[str, object]]:
    staging_dir = atlas_root / "spectre"
    ras_orientation = axcodes2ornt(("R", "A", "S"))
    staged_files = sorted(staging_dir.glob("*.nii.gz"))
    if not staged_files:
        raise RuntimeError(f"No staged NIfTI files found in {staging_dir}")

    inventory: list[dict[str, object]] = []
    destinations: set[Path] = set()
    for staged in staged_files:
        stem = staged.name.removesuffix(".nii.gz")
        sources = [
            path
            for pattern in (f"{stem}.raw", f"{stem}.img")
            for path in atlas_root.rglob(pattern)
        ]
        if len(sources) != 1:
            raise RuntimeError(f"Expected one legacy image for {staged.name}, found {sources}")
        destination = sources[0].parent / staged.name
        if (destination.exists() and destination != staged) or destination in destinations:
            raise FileExistsError(f"Destination already exists or is duplicated: {destination}")
        destinations.add(destination)

        image = nib.load(str(staged))
        original_orientation = io_orientation(image.affine)
        if np.isnan(original_orientation).any():
            raise ValueError(f"Indeterminate orientation in {staged}")
        transform = ornt_transform(original_orientation, ras_orientation)
        canonical = nib.as_closest_canonical(image)
        source_values = np.asanyarray(image.dataobj)
        expected_values = apply_orientation(source_values, transform)
        canonical_values = np.asanyarray(canonical.dataobj)
        if not np.array_equal(expected_values, canonical_values):
            raise ValueError(f"Voxel values changed while normalizing {staged}")
        expected_affine = image.affine @ inv_ornt_aff(transform, image.shape)
        if not np.allclose(canonical.affine, expected_affine, rtol=0, atol=1e-5):
            raise ValueError(f"Physical mapping changed while normalizing {staged}")
        if nib.aff2axcodes(canonical.affine) != ("R", "A", "S"):
            raise ValueError(f"Failed to make {staged} RAS-oriented")

        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(
            f".{destination.name.removesuffix('.nii.gz')}.tmp.nii.gz"
        )
        nib.save(canonical, str(temporary))
        temporary.replace(destination)
        inventory.append(
            {
                "staged_filename": staged.name,
                "source_payload": str(sources[0].relative_to(atlas_root)),
                "source_payload_sha256": sha256(sources[0]),
                "mipav_xml": (
                    str(sources[0].with_suffix(".xml").relative_to(atlas_root))
                    if sources[0].with_suffix(".xml").exists()
                    else None
                ),
                "mipav_xml_sha256": (
                    sha256(sources[0].with_suffix(".xml"))
                    if sources[0].with_suffix(".xml").exists()
                    else None
                ),
                "nifti_input_sha256": sha256(staged),
                "nifti_output": str(destination.relative_to(atlas_root)),
                "nifti_output_sha256": sha256(destination),
                "input_shape": list(image.shape),
                "output_shape": list(canonical.shape),
                "input_dtype": str(image.get_data_dtype()),
                "output_dtype": str(canonical.get_data_dtype()),
                "input_axis_codes": list(nib.aff2axcodes(image.affine)),
                "output_axis_codes": list(nib.aff2axcodes(canonical.affine)),
                "input_affine_ras": _affine_list(image.affine),
                "output_affine_ras": _affine_list(canonical.affine),
                "voxel_values_exact_after_axis_reorientation": True,
                "physical_mapping_preserved_from_supplied_nifti": True,
            }
        )

    # Update descriptors after all image outputs have passed orientation,
    # value, geometry, and checksum checks.
    for descriptor in atlas_root.rglob("*.txt"):
        original_text = descriptor.read_text(encoding="utf-8")
        updated_text = original_text.replace(".raw", ".nii.gz").replace(".img", ".nii.gz")
        if updated_text != original_text:
            descriptor.write_text(updated_text, encoding="utf-8")

    audit_path = atlas_root / "CONVERSION_AUDIT.json"
    audit_path.write_text(
        json.dumps(
            {
                "description": "User-supplied NIfTI conversions normalized to RAS.",
                "method": "NiBabel axis permutation/flip only; no interpolation.",
                "images": inventory,
                "known_source_metadata_caveats": [
                    "Some source MIPAV XML files declare Unknown subject-axis orientation; their supplied NIfTI affine was used and preserved.",
                    "OAS1 rotated raw volumes have a Little endian XML declaration, while big-endian decoding exactly matches the supplied NIfTI values; preserve this discrepancy for source review.",
                    "Analyze-derived OAS1 NIfTI affines do not establish a trusted physical origin; verify their intended origin before scientific use.",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return inventory


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1] / "data" / "atlases"
    rows = normalize_atlases(root)
    print(f"Normalized {len(rows)} atlas NIfTIs to RAS and updated descriptors.")
