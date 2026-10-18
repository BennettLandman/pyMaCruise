"""Verify installed asset completeness without accessing legacy checkouts.

Run ``python -m pymacruise.asset_inventory`` after installation. Scientific
equivalence and native-atlas coordinate conversion are separate validation.
"""
from __future__ import annotations

import hashlib
from importlib.resources import files
import json
from pathlib import Path

import nibabel as nib

from pymacruise.atlas_resources import (
    get_bundled_atlas_descriptor, get_bundled_atlas_root, read_bundled_atlas_source,
)


def _check_manifest(root: Path) -> int:
    entries = {}
    for line in (root / "SHA256SUMS").read_text().splitlines():
        expected, relative = line.split("  ", 1)
        unresolved = root / relative
        path = unresolved.resolve()
        if not path.is_relative_to(root.resolve()) or unresolved.is_symlink() or not path.is_file():
            raise RuntimeError(f"Missing or unsafe bundled asset: {relative}")
        if relative in entries:
            raise RuntimeError(f"Duplicate bundled asset checksum: {relative}")
        entries[relative] = expected
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Bundled asset checksum mismatch: {path}")
    actual = {p.relative_to(root).as_posix() for p in root.rglob("*")
              if p.is_file() and p.name != "SHA256SUMS"}
    if actual != set(entries):
        raise RuntimeError(f"Asset manifest file set differs: {actual.symmetric_difference(entries)}")
    return len(entries)


def verify_bundled_assets() -> dict:
    """Check every stored file and decode native sources and runtime tables."""
    root = get_bundled_atlas_root().resolve()
    atlas_files = _check_manifest(root)
    runtime_root = files("pymacruise").joinpath("data")
    # A wheel exposes paths; source checkouts also use a real package tree.
    runtime_root = Path(str(runtime_root))
    table_entries = (runtime_root / "SHA256SUMS").read_text().splitlines()
    table_names = [line.split("  ", 1)[1] for line in table_entries]
    if len(table_names) != 2 or set(table_names) != {
        "critical186LUT.raw.gz", "jist_isosurface_connectivity18_cases.json"
    }:
        raise RuntimeError("Runtime table manifest must cover both required tables exactly once")
    for line in table_entries:
        expected, relative = line.split("  ", 1)
        path = (runtime_root / relative).resolve()
        if not path.is_relative_to(runtime_root.resolve()) or not path.is_file():
            raise RuntimeError(f"Missing or unsafe runtime table: {relative}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Runtime table checksum mismatch: {path}")
    from pymacruise.toads_critical_lut import load_critical186_lut
    from pymacruise.surfaces import _jist_connectivity18_cases
    assert load_critical186_lut().shape == (1 << 26,)
    assert len(_jist_connectivity18_cases()) == 256

    nifti_files = list(root.rglob("*.nii.gz"))
    for path in nifti_files:
        image = nib.load(path)
        if len(image.shape) != 3 or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise RuntimeError(f"Expected a 3D RAS atlas image: {path}")
    descriptor_count = 0
    for path in root.rglob("*.txt"):
        for line in path.read_text().splitlines():
            if line.startswith("File:"):
                relative = line.partition(":")[2].strip()
                payload = (path.parent / relative).resolve()
                if not payload.is_relative_to(root) or not payload.is_file():
                    raise RuntimeError(f"Missing descriptor payload {path}: {relative}")
        descriptor_count += 1
    for descriptor in ("cruise-atlas-10obj-toads2008.txt",
                       "lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"):
        get_bundled_atlas_descriptor(descriptor)

    native = json.loads((root / "ibsr11objV3/SOURCE_MANIFEST.json").read_text())
    for filename in native["files"]:
        read_bundled_atlas_source("ibsr11objV3", filename)
    for reference in native["referenced_payloads"]:
        if Path(reference).name not in native["files"]:
            raise RuntimeError(f"Missing native PLUS descriptor payload: {reference}")
    return {
        "atlas_files_checked": atlas_files, "ras_nifti_files_checked": len(nifti_files),
        "runtime_tables_checked": len(table_entries), "descriptors_checked": descriptor_count,
        "ibsr11objV3_source_files_checked": len(native["files"]),
        "ibsr11objV3_referenced_payloads_checked": len(native["referenced_payloads"]),
        "ibsr11objV3_runtime_status": native["runtime_status"],
    }


if __name__ == "__main__":
    print(json.dumps(verify_bundled_assets(), indent=2, sort_keys=True))
