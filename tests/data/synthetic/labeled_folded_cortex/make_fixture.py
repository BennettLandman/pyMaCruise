"""Generate a RAS folded cortex with explicit MaSeg labels and TGDM inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
SHAPE = (65, 65, 65)
CENTER = np.asarray([32.0, 32.0, 32.0], dtype=np.float32)
AFFINE = np.eye(4, dtype=np.float64)


def fixture_arrays() -> tuple[dict[str, np.ndarray], dict[str, object]]:
    x, y, z = np.indices(SHAPE, dtype=np.float32)
    dx, dy, dz = x - CENTER[0], y - CENTER[1], z - CENTER[2]
    radius = np.sqrt(dx * dx + dy * dy + dz * dz)
    theta = np.arctan2(dy, dx)
    phi = np.arctan2(np.sqrt(dx * dx + dy * dy), dz)
    # Low-order lobes create cortex-like gyri/sulci without self-intersection.
    fold = 1.8 * np.sin(phi) ** 2 * np.cos(6.0 * theta)
    fold += 0.9 * np.sin(phi) ** 2 * np.cos(4.0 * phi + 0.35)
    fold += 0.35 * np.sin(phi) ** 4 * np.sin(3.0 * theta + 0.2)
    outer_radius = np.float32(22.0) + fold
    inner_radius = outer_radius - np.float32(5.0)
    brain_radius = outer_radius + np.float32(2.0)

    wm = radius <= inner_radius
    gm = (radius > inner_radius) & (radius <= outer_radius)
    csf = (radius > outer_radius) & (radius <= brain_radius)
    brain = radius <= brain_radius
    labels = np.zeros(SHAPE, dtype=np.uint8)
    labels[csf] = 5
    labels[wm] = 44
    labels[gm] = 100

    t1 = np.zeros(SHAPE, dtype=np.float32)
    t1[wm] = 110.0
    t1[gm] = 145.0
    t1[csf] = 35.0
    membership = np.where(brain, 255, 0).astype(np.uint8)
    brain_bound = np.where(brain, -1.0, 1.0).astype(np.float32)
    white_matter = np.where(wm, 255.0, 0.0).astype(np.float32)
    gray_matter = np.where(gm, 255.0, 0.0).astype(np.float32)
    # Hold the optional skeleton force off in this TGDM fixture. MaACE is
    # validated separately; this case isolates the three source surfaces.
    ace_skeleton = np.zeros(SHAPE, dtype=np.uint8)
    truth_inner = radius - inner_radius
    truth_outer = radius - outer_radius
    arrays = {
        "t1.nii.gz": t1,
        "brain_mask.nii.gz": membership,
        "tgdm_brain_bound.nii.gz": brain_bound,
        "maseg_labels.nii.gz": labels,
        "initial_membership.nii.gz": membership,
        "white_matter.nii.gz": white_matter,
        "gray_matter.nii.gz": gray_matter,
        "ace_gray_matter.nii.gz": gray_matter.copy(),
        "ace_skeleton.nii.gz": ace_skeleton,
        "ground_truth_inner_levelset.nii.gz": truth_inner.astype(np.float32),
        "ground_truth_outer_levelset.nii.gz": truth_outer.astype(np.float32),
    }
    spec: dict[str, object] = {
        "description": "Deterministic folded cortical ribbon with explicit MaSeg tissue labels; entirely synthetic.",
        "shape_voxels": list(SHAPE),
        "center_voxel_ras": CENTER.tolist(),
        "affine_ras_mm": AFFINE.tolist(),
        "orientation": "RAS",
        "voxel_size_mm": [1.0, 1.0, 1.0],
        "outer_radius_mm": "22 + 1.8*sin(phi)^2*cos(6*theta) + 0.9*sin(phi)^2*cos(4*phi+0.35) + 0.35*sin(phi)^4*sin(3*theta+0.2)",
        "ground_truth_cortical_ribbon_thickness_mm": 5.0,
        "label_values": {"background": 0, "sulcal_csf": 5, "white_matter": 44, "cortical_gray_matter": 100},
        "image_intensities": {"background": 0, "sulcal_csf": 35, "white_matter": 110, "cortical_gray_matter": 145},
        "supports": {
            "initial_membership": "255 within the binary brain mask; 0 outside",
            "tgdm_brain_bound": "inside-negative/outside-positive unit-sign field derived from brain_mask for GenericTGDM bound semantics",
            "white_matter": "255 inside the analytic WM boundary; 0 elsewhere",
            "gray_matter": "255 in the analytic cortical ribbon; 0 elsewhere",
            "ace_gray_matter": "same support and values as gray_matter",
            "ace_skeleton": "zero everywhere to isolate nested-surface forces; MaACE is covered by separate source comparisons",
            "ground_truth_inner_levelset": "radius minus analytic inner boundary; negative inside WM",
            "ground_truth_outer_levelset": "radius minus analytic outer boundary; negative inside cortex",
        },
        "selected_maCRUISE_surface_evolution": {
            "initial_iso_level": 126.9,
            "white_matter_iso": 0.6,
            "gray_matter_iso": 0.4,
            "initial_smoothing_curvature": 1.0,
            "initial_smoothing_iterations": 6,
            "inner_curvature": 0.2,
            "inner_pressure": 1.0,
            "inner_iterations": 4,
            "central_curvature": 0.15,
            "central_pressure": 1.5,
            "central_external": 1.0,
            "central_iterations": 7,
            "central_shrink_iterations": 1,
            "outer_curvature": 1.0,
            "outer_pressure": 1.0,
            "outer_iterations": 7,
            "provenance": "TGDMBaselineParameters selected by run_ToadsCruise; WM iso is overridden to 0.6 in the base MaCRUISE core.",
        },
        "ribbon_quality_checks": [
            "inner, central, and outer meshes are finite, closed, manifold, one component, Euler characteristic 2, and RAS-mm",
            "each adjacent pair is strictly nested with no zero/negative sampled ribbon thickness",
            "WM label 44 is inside the inner surface, GM label 100 lies between inner and outer, and CSF label 5 lies outside the outer surface, with voxel-boundary tolerance stated in the report",
            "Java/Python final inner/central/outer surfaces meet Dice >= 0.97, symmetric mean distance <= 0.01 mm, and Hausdorff <= 1 mm",
            "positive finite thickness is produced wherever inner/outer surfaces have valid opposing closest points",
        ],
    }
    return arrays, spec


def main() -> None:
    arrays, spec = fixture_arrays()
    for name, values in arrays.items():
        image = nib.Nifti1Image(values, AFFINE)
        image.set_qform(AFFINE, code=1)
        image.set_sform(AFFINE, code=1)
        nib.save(image, ROOT / name)
    (ROOT / "spec.json").write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    payloads = sorted(ROOT.glob("*.nii.gz")) + [ROOT / "spec.json"]
    (ROOT / "SHA256SUMS").write_text(
        "".join(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n" for path in payloads),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
