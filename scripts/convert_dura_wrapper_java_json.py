"""Convert selected Dura wrapper Java JSON arrays into RAS NIfTI fixtures."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURE = ROOT / "tests/data/synthetic/dura_attached_no_dura_41"
AFFINE = np.eye(4, dtype=np.float64)
FLOAT_ARRAYS = (
    "t1", "wm", "gm", "wm_tc", "union", "gm_tc", "inner_distance",
    "outer_distance", "level_set", "stripped_t1", "corrected_wm", "corrected_gm",
    "brain_mask_level_set",
)


def save(path: Path, values: np.ndarray, *, byte: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dtype = np.uint8 if byte else np.float32
    nib.save(nib.Nifti1Image(np.asarray(values, dtype=dtype), AFFINE), path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_json", type=Path)
    parser.add_argument("--fixture-dir", type=Path, default=DEFAULT_FIXTURE)
    args = parser.parse_args()
    fixture = args.fixture_dir
    reference = json.loads(args.reference_json.read_text())
    dimensions = reference.get("dimensions_xyz")
    if not isinstance(dimensions, list) or len(dimensions) != 3 or len(set(dimensions)) != 1:
        raise ValueError("expected selected Java Dura wrapper output for a cubic grid")
    shape = tuple(int(value) for value in dimensions)
    source_hashes = {
        "MedicAlgorithmRemoveDura09": reference["medic_algorithm_remove_dura09_sha256"],
        "AlgorithmTopologyCorrection": reference["algorithm_topology_correction_sha256"],
        "GenericTGDM": reference["generic_tgdm_sha256"],
        "DynamicPressureTGDM": reference["dynamic_pressure_tgdm_sha256"],
        "DynamicPressureCortexProfile": reference["cortex_profile_sha256"],
        "DistanceField": reference["distance_field_sha256"],
    }
    if source_hashes.get("MedicAlgorithmRemoveDura09") != "1ed180daf4a4334c995ad310828ec91cfd5593413eec17ab15404a421d7805a4":
        raise ValueError("unexpected MedicAlgorithmRemoveDura09 source revision")

    for case in reference["cases"]:
        case_root = fixture / case["name"]
        inputs = case_root / "inputs"
        expected = case_root / "java_reference"
        arrays = case["arrays"]
        for name in ("t1", "wm", "gm", "paint"):
            save(inputs / f"{name}.nii.gz", np.asarray(arrays[name]).reshape(shape), byte=(name == "paint"))
        for name in FLOAT_ARRAYS:
            save(expected / f"{name}.nii.gz", np.asarray(arrays[name]).reshape(shape))
        save(expected / "wm_topology.nii.gz", np.asarray(arrays["wm_tc"]).reshape(shape))
        save(expected / "gm_topology.nii.gz", np.asarray(arrays["gm_tc"]).reshape(shape))
        save(expected / "brain_mask.nii.gz", np.asarray(arrays["brain_mask"]).reshape(shape), byte=True)
        for index, values in enumerate(case["trace_levelsets"], start=1):
            save(expected / f"tgdm_pass_{index:02d}.nii.gz", np.asarray(values[f"pass_{index-1}"]).reshape(shape))
        trace = {
            "trace_steps": case["trace_steps"],
            "trace_levelset_sha256": case["trace_levelset_sha256"],
            "source_sha256": source_hashes,
            "case": case["name"],
            "dimensions_xyz": reference["dimensions_xyz"],
        }
        expected.mkdir(parents=True, exist_ok=True)
        (expected / "trace.json").write_text(json.dumps(trace, indent=2) + "\n")

    hashes = {
        str(path.relative_to(fixture)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(fixture.rglob("*.nii.gz"))
    }
    (fixture / "SHA256SUMS").write_text("".join(f"{digest}  {path}\n" for path, digest in hashes.items()))
    provenance = {
        "scope": "synthetic selected-wrapper-sequence comparison",
        "provenance": "analytic scaled synthetic attached-dura and no-dura cases; no subject data",
        "dimensions_xyz": dimensions,
        "scale": reference["scale"],
        "affine_ras_mm": AFFINE.tolist(),
        "source_sha256": source_hashes,
        "harness_sha256": hashlib.sha256((ROOT / "tests/java_reference/DuraAttachedNoDuraReferenceHarness.java").read_bytes()).hexdigest(),
        "cases": [case["name"] for case in reference["cases"]],
        "parameters": {
            "membership_level": reference["membership_level"],
            "mask_translation": reference["mask_translation"],
            "level_set_translation": reference["level_set_translation"],
            "wm_topology_regularization": 0.05,
            "gm_topology_regularization": float(reference.get("gm_regularization", 0.02)),
            "tgdm_passes": 5,
        },
        "java_nifti_sha256": hashes,
        "reference_scope_note": "Java harness follows MedicAlgorithmRemoveDura09.execute() using selected topology, distance, and TGDM classes; it is a wrapper-sequence reference, not a MIPAV/JIST application run.",
    }
    (fixture / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (fixture / "fixture_hashes.json").write_text(json.dumps(hashes, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"fixture": str(fixture), "cases": provenance["cases"], "nifti_count": len(hashes)}, indent=2))


if __name__ == "__main__":
    main()
