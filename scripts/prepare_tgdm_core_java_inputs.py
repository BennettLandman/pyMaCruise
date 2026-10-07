#!/usr/bin/env python3
"""Capture the synthetic base-run TGDM inputs as a reproducible Java fixture."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TARGET = REPO_ROOT / "tests/data/synthetic/tgdm_core_brain_resolution/java_inputs"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("core_dir", type=Path, help="completed synthetic toads-cruise-core output directory")
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    args = parser.parse_args()
    args.target.mkdir(parents=True, exist_ok=True)

    mapping = {
        "initial_membership": "wmfill_topology_corrected.nii.gz",
        "white_matter": "wmfill_cropped.nii.gz",
        "gray_matter": "gm_cropped.nii.gz",
        "ace_gray_matter": "ace_enhanced_gm.nii.gz",
        "ace_skeleton": "ace_thinned_skeleton.nii.gz",
        "gvf": "gvf.nii.gz",
    }
    images: dict[str, nib.spatialimages.SpatialImage] = {}
    for name, filename in mapping.items():
        image = nib.load(str(args.core_dir / filename))
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"{filename} is not RAS oriented")
        if name == "gvf":
            if image.ndim != 4 or image.shape[-1] != 3:
                raise ValueError("GVF must be a component-last 3D vector image")
        elif image.ndim != 3:
            raise ValueError(f"{filename} must be scalar 3D")
        images[name] = image

    reference = images["initial_membership"]
    for name, image in images.items():
        if image.shape[:3] != reference.shape[:3] or not np.allclose(image.affine, reference.affine):
            raise ValueError(f"{name} is not aligned to the initial membership grid")
        target = args.target / f"{name}.nii.gz"
        data = np.asarray(image.dataobj)
        nib.save(nib.Nifti1Image(data, image.affine, image.header.copy()), str(target))

    source_provenance = args.core_dir / "provenance.json"
    provenance = json.loads(source_provenance.read_text()) if source_provenance.exists() else {}
    manifest = {
        "synthetic_only": True,
        "source": "repository-generated synthetic atlas-derived base workflow; no human scan",
        "source_stage": str(args.core_dir),
        "orientation": "RAS",
        "shape": list(reference.shape),
        "affine": reference.affine.tolist(),
        "inputs": {name: {"file": f"{name}.nii.gz", "sha256": sha256(args.target / f"{name}.nii.gz")} for name in mapping},
        "selected_wrapper_parameters": provenance.get("parameters", {}).get("tgdm", {}),
    }
    (args.target / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Prepared {len(mapping)} synthetic RAS TGDM inputs in {args.target}")


if __name__ == "__main__":
    main()
