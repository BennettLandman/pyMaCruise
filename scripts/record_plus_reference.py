#!/usr/bin/env python3
"""Inventory a completed original PLUS driver run for the product comparator.

Reads only original outputs and the shared testing kit. Does not read or use
Python segmentation results. Run only after the original process exits zero.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re

import nibabel as nib
import numpy as np
from compare_external_plus_run import PRODUCTS, VOLUMES, _image, _regular, sha256


def payload_hashes(descriptor: Path) -> dict[str, str]:
    result = {}
    for value in re.findall(r"^File:\s*(.+?)\s*$", descriptor.read_text(), re.MULTILINE):
        path = descriptor.parent / value
        key = path.name
        if key in result:
            # Some descriptor structures intentionally share a payload.
            if result[key] != sha256(_regular(path, "atlas payload")):
                raise ValueError(f"Ambiguous atlas payload basename: {key}")
        result[key] = sha256(_regular(path, "atlas payload"))
    if not result:
        raise ValueError(f"No payload references in {descriptor}")
    return result


def record(args):
    if args.exit_code != 0:
        raise ValueError("Keep failed logs; a failed original run cannot be recorded as complete")
    if any(not value.strip() or 'REPLACE' in value for value in (args.runtime, args.source_revision, args.command)):
        raise ValueError("Supply actual runtime, source revision and launch command, replacing template text")
    root = args.source_root.resolve()
    kit = json.loads(args.kit.read_text())
    outputs = json.loads((root / "driver_outputs.json").read_text())
    input_dir = Path(outputs["input_dir"]).resolve()
    input_hashes = {}
    for name, expected in kit["inputs"].items():
        path = Path(expected["path"])
        if input_dir != path.parent or sha256(path) != expected["sha256"]:
            raise ValueError(f"Original input directory/hash differs from kit: {name}")
        input_hashes[name] = expected["sha256"]
    descriptor = _regular(Path(outputs["atlasem_descriptor"]), "original AtlasEM descriptor")
    if descriptor != Path(kit["source_atlasem_descriptor"]).resolve():
        raise ValueError("Original driver did not use the kit's selected native AtlasEM descriptor")
    if sha256(descriptor) != kit["source_atlasem_descriptor_sha256"]:
        raise ValueError("Native AtlasEM descriptor changed since kit preparation")
    products = {name: Path(outputs["refined_segmentation" if name == "segmentation" else name]).resolve()
                for name in PRODUCTS}
    for name, path in products.items():
        _regular(path, name)
        if not path.is_relative_to(root):
            raise ValueError(f"Original product is outside the returned run: {name}")
        if name in VOLUMES:
            _image(path)
    labels = np.unique(np.asarray(nib.load(kit["inputs"]["maseg"]["path"]).dataobj))
    if not np.array_equal(labels, np.rint(labels)):
        raise ValueError("Starting MaSeg is not integer-labelled")
    lesion = _regular(args.lesion_descriptor.resolve(), "original LesionTOADS descriptor")
    input_hashes["atlasem_descriptor"] = sha256(descriptor)
    input_hashes["atlas_descriptor"] = sha256(lesion)
    log = _regular(args.log.resolve(), "original launch log")
    if log.stat().st_size == 0:
        raise ValueError("Original log is empty")
    info = {
        "status": "complete", "exit_code": args.exit_code,
        "runtime": args.runtime, "source_revision": args.source_revision, "command": args.command,
        "atlasem_enabled": True, "atlasem_descriptor_sha256": sha256(descriptor),
        "atlasem_asset_sha256": payload_hashes(descriptor),
        "lesion_atlas_asset_sha256": payload_hashes(lesion),
        "input_sha256": input_hashes,
        "products_scope": "independently_generated_selected_plus_products",
        "preprocessing_boundary": "shared raw T1/FLAIR/MaSeg and user brain mask; original Aladin/N4/AtlasEM/Dura/LesionTOADS scientific stages; no cross-pipeline scientific intermediates",
        "surface_frame": "T1 voxel indices" if args.mesh_frame == "voxel" else "RAS+ mm",
        "required_segmentation_labels": [int(value) for value in labels if value != 0],
        "active_product_attempts": {name: args.attempt for name in PRODUCTS},
        "product_sha256": {name: sha256(path) for name, path in products.items()},
        "source_log": {"path": str(log), "sha256": sha256(log)},
        "testing_kit": str(args.kit.resolve()),
        "limitations": ["Runtime, revision, launch command and mesh frame are supplied by the operator and require review.",
                        "Atlas byte differences are expected for native Java versus converted Python representations; numerical atlas equivalence remains to be checked."],
    }
    (root / "product-map.json").write_text(json.dumps({name: str(path.relative_to(root)) for name, path in products.items()}, indent=2) + "\n")
    (root / "reference-info.json").write_text(json.dumps(info, indent=2) + "\n")
    return root / "reference-info.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("kit", "source-root", "lesion-descriptor", "log"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("runtime", "source-revision", "command"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--exit-code", type=int, required=True)
    parser.add_argument("--mesh-frame", choices=("voxel", "ras-mm"), required=True)
    parser.add_argument("--attempt", default="original-attempt-01")
    args = parser.parse_args()
    try:
        print(record(args))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"Cannot record original PLUS run: {exc}\n")


if __name__ == "__main__":
    main()
