from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "tests"))
from test_lesion_post_thinning_gm_membership import SPEC, _compose  # noqa: E402
from pymacruise.lesion_post_thinning_membership import (  # noqa: E402
    run_post_thinning_standard13_membership_update,
)


ROOT = Path(__file__).resolve().parent
SOURCE_FIXTURE = REPO / "tests/data/synthetic/lesion_initial_thinning_composition/initial_thinning_inputs_ras.nii.gz"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    packed, state, _, thinning = _compose()
    update = run_post_thinning_standard13_membership_update(
        thinning,
        state["images"],
        SPEC["class_centroids"],
        state["lesion"],
        state["black"],
        state["outlier"],
        SPEC["lesion_centroid"],
        (0.0, 0.0),
        state["shape_values"],
        SPEC["shape_priors"],
        SPEC["smoothing_priors"],
        state["direct_relations"],
        packed[..., 63],
        state["ventricle_distances"],
        state["gray_matter_distances"],
        state["brainstem_distances"],
        SPEC["object_types"],
        centroid_options={
            "mode": "none",
            "modalities": ("T1_MPRAGE", "FLAIR"),
            "lesion_power_exponents": (0.4, 1.2),
        },
        channel_weight_options={"outlier_distance": SPEC["outlier_energy"]},
        outlier_energy=SPEC["outlier_energy"],
        lesion_power_exponent=SPEC["lesion_power_exponent"],
        flair_channel=1,
    )
    membership = update.membership.class_memberships
    if not np.isfinite(membership).all():
        raise ValueError("selected registration requires finite membership maps")
    shape_priors = np.asarray(state["shape_values"], dtype=np.float32)
    combined = np.concatenate((membership, shape_priors), axis=-1).astype(np.float32)
    source_image = nib.load(str(SOURCE_FIXTURE))
    output = ROOT / "registration_inputs_ras.nii.gz"
    nib.save(nib.Nifti1Image(combined, source_image.affine), str(output))
    metadata = {
        "description": "P39am GMflag=true class memberships and fixed class shape priors from the folded synthetic P39ak RAS fixture",
        "dimensions": list(combined.shape),
        "frame_order": [f"membership_{i}" for i in range(13)]
        + [f"shape_prior_{i}" for i in range(13)],
        "affine_ras_mm": source_image.affine.tolist(),
        "initial_transform": [0.03, -0.02, 0.01, 0.1, -0.1, 0.05],
        "image_resolutions_mm": [1.0, 1.0, 1.0],
        "template_resolutions_mm": [1.0, 1.0, 1.0],
        "subsample": 2,
        "init_align_iter": 10,
        "source_call_order": ["computeMemberships(GMflag=true)", "energyRatio", "registerShapes"],
        "membership_nan_objective_does_not_gate_registration": True,
        "source_fixture_sha256": sha256(SOURCE_FIXTURE),
        "fixture_sha256": sha256(output),
        "membership_count": int(membership.size),
        "membership_nan_count": int(np.count_nonzero(np.isnan(membership))),
        "shape_support_per_class": [int(np.count_nonzero(shape_priors[..., c] > 0)) for c in range(13)],
    }
    (ROOT / "fixture.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
