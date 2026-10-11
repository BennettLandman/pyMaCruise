#!/usr/bin/env python3
"""Compare preserved CombineMAandTOADS Java source with the plus Python stage."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy import ndimage


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    root = args.reference_dir
    provenance = json.loads((root / "topology-provenance.json").read_text())
    shape = tuple(provenance["shape_xyz"])
    rows = {}
    for name in ("gm", "wm", "mask"):
        with gzip.open(root / f"java-{name}.raw.gz", "rb") as stream:
            header = tuple(int.from_bytes(stream.read(4), "big") for _ in range(3))
            java = np.frombuffer(stream.read(), dtype=">f4").astype(np.float32)
        if header != shape or java.size != int(np.prod(shape)):
            raise SystemExit(f"{name}: Java output shape mismatch")
        java = java.reshape(shape, order="C")
        python = np.load(root / f"python-tissue-{name}.npy").astype(np.float32)
        diff = np.abs(java.astype(np.float64) - python.astype(np.float64))
        rows[name] = {
            "exact": bool(np.array_equal(java, python)),
            "changed_voxels": int(np.count_nonzero(java != python)),
            "voxel_count": int(java.size),
            "maximum_absolute_error": float(diff.max(initial=0)),
            "mean_absolute_error": float(diff.mean()),
        }

    tissue = {
        name: np.load(root / f"python-input-{name}.npy")
        for name in (
            "atlas-gm", "atlas-wm", "brain-mask", "toads-gm",
            "toads-wm-fill", "toads-wm-mask",
        )
    }
    maseg = np.load(root / "python-input-maseg.npy")
    label_100 = maseg == 100
    wm_int = np.trunc(tissue["atlas-wm"]).astype(np.int32)
    eroded_wm = np.zeros(wm_int.shape, dtype=np.int32)
    if all(length >= 3 for length in wm_int.shape):
        six_neighbor = ndimage.generate_binary_structure(3, 1)
        min_neighbors = ndimage.minimum_filter(wm_int, footprint=six_neighbor, mode="constant", cval=0)
        interior = tuple(slice(1, -1) for _ in range(3))
        eroded_wm[interior] = min_neighbors[interior]
    inside = np.trunc(tissue["brain-mask"]).astype(np.int32) > 0
    confident_wm = eroded_wm > 0
    confident_gm = (
        (tissue["toads-gm"] > 0)
        & (tissue["toads-gm"] < 0.255)
        & (np.trunc(tissue["atlas-gm"]).astype(np.int32) > 0)
    )
    branch_masks = {
        "outside_mask": ~inside,
        "confident_wm": inside & confident_wm,
        "confident_gm": inside & ~confident_wm & confident_gm,
        "fallback": inside & ~confident_wm & ~confident_gm,
    }
    label_100_diagnostics = {
        "voxel_count": int(np.count_nonzero(label_100)),
        "fusion_mask_support_voxels": int(np.count_nonzero(label_100 & inside)),
        "atlas_gm_support_voxels": int(np.count_nonzero(label_100 & (tissue["atlas-gm"] > 0))),
        "toads_gm_nonzero_voxels": int(np.count_nonzero(label_100 & (tissue["toads-gm"] != 0))),
        "toads_wm_fill_nonzero_voxels": int(np.count_nonzero(label_100 & (tissue["toads-wm-fill"] != 0))),
        "branch_counts": {
            name: int(np.count_nonzero(label_100 & mask)) for name, mask in branch_masks.items()
        },
        "fused_gm_nonzero_voxels": int(
            np.count_nonzero(label_100 & (np.load(root / "python-tissue-gm.npy") != 0))
        ),
        "fused_wm_nonzero_voxels": int(
            np.count_nonzero(label_100 & (np.load(root / "python-tissue-wm.npy") != 0))
        ),
    }
    report = {
        "synthetic_only": True,
        "comparison": "preserved selected Java CombineMAandTOADS.execute() vs Python combine_multi_atlas_toads on workflow crop inputs",
        "scope": "held-fixed six atlas/TOADS input images from the saved plus run; selected default ErodeMethod=3D and GMconfident=0.001",
        "workflow_run_signature": provenance["workflow_run_signature"],
        "shape_xyz": list(shape),
        "coordinate_system": "RAS",
        "brain_mask_override_source": provenance.get("brain_mask_override_source"),
        "label_100_diagnostics": label_100_diagnostics,
        "outputs": rows,
        "all_outputs_exact": all(record["exact"] for record in rows.values()),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))
    if not report["all_outputs_exact"]:
        raise SystemExit("Java and Python plus tissue products differ")


if __name__ == "__main__":
    main()
