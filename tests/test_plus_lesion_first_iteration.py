from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.plus_lesion_atlas_init import run_plus_lesion_initial_atlas_alignment_checkpoint
from pymacruise.plus_lesion_crop import run_plus_lesion_crop_checkpoint
import pymacruise.plus_lesion_first_iteration as first_stage
from pymacruise.plus_lesion_registered_recrop import run_plus_lesion_registered_recrop_checkpoint
from pymacruise.lesion_registered_first_iteration import (
    run_lesion_toads_registered_first_iteration,
)
from pymacruise.lesion_registered_first_iteration_branch import (
    run_lesion_toads_registered_first_iteration_branch,
)
from pymacruise.lesion_atlas_priors import (
    load_lesion_atlas_shape_priors,
    load_lesion_atlas_topology,
)
from pymacruise.lesion_standard13_voxel import STANDARD13_TYPES


ROOT = Path(__file__).parent
SYNTHETIC = ROOT / "data/synthetic/plus_lesion_atlas_init"
ATLAS = SYNTHETIC / "mini-lesion-toads-atlas.txt"


def _run_recrop(tmp_path: Path):
    crop = run_plus_lesion_crop_checkpoint(
        SYNTHETIC / "synthetic_t1_flair.nii.gz",
        SYNTHETIC / "synthetic_topology.nii.gz",
        tmp_path / "crop",
        padding=1,
    )
    alignment = run_plus_lesion_initial_atlas_alignment_checkpoint(
        crop.cropped_path,
        crop.manifest_path,
        ATLAS,
        tmp_path / "alignment",
        init_align_iterations=0,
        pyramid_levels=2,
    )
    return run_plus_lesion_registered_recrop_checkpoint(
        crop.normalized_path,
        crop.manifest_path,
        alignment.manifest_path,
        ATLAS,
        tmp_path / "registered-recrops",
    )


def test_first_iteration_checkpoint_saves_source_ordered_membership_and_thinning(
    tmp_path: Path, monkeypatch
) -> None:
    recrop = _run_recrop(tmp_path)
    with pytest.raises(ValueError, match="initial membership objective must be finite"):
        first_stage.run_plus_lesion_first_iteration_checkpoint(
            recrop.manifest_path,
            ATLAS,
            tmp_path / "source-faithful",
            iterations=1,
            first_thinning_limit=1000.0,
        )
    failed = json.loads((tmp_path / "source-faithful/lesion-first-iteration.json").read_text())
    assert failed["status"] == "failed"
    assert failed["failure"]["source_zero_variance_nan"] is True
    assert failed["failure"]["initial_objective_finite"] is False
    assert failed["failure"]["nonfinite_initial_membership_values"] == 0
    assert failed["failure"]["nonfinite_initial_energy_values"] > 0
    assert failed["failure"]["zero_lesion_standard_deviation_channels"] == 2

    captured = {}
    source_initial = first_stage.run_lesion_toads_registered_initial_membership

    def capture_initial(*args, **kwargs):
        captured["initial"] = source_initial(*args, **kwargs)
        return captured["initial"]

    monkeypatch.setattr(
        first_stage, "run_lesion_toads_registered_initial_membership", capture_initial
    )
    result = first_stage.run_plus_lesion_first_iteration_checkpoint(
        recrop.manifest_path,
        ATLAS,
        tmp_path / "first-iteration",
        iterations=1,
        first_thinning_limit=1000.0,
        source_zero_variance_nan=False,
    )
    modalities = nib.load(str(recrop.modalities_path))
    memberships = nib.load(str(result.memberships_path))
    segmentation = nib.load(str(result.segmentation_path))
    skeleton = nib.load(str(result.skeleton_path))
    availability = nib.load(str(result.availability_path))
    ordering = nib.load(str(result.ordering_path))
    assert memberships.shape == (*modalities.shape[:3], 17)
    assert segmentation.shape == modalities.shape[:3]
    assert skeleton.shape == modalities.shape[:3]
    assert availability.shape == (*modalities.shape[:3], 13)
    assert ordering.shape == modalities.shape[:3]
    assert all(
        nib.aff2axcodes(image.affine) == ("R", "A", "S")
        for image in (memberships, segmentation, skeleton, availability, ordering)
    )
    assert np.isfinite(np.asarray(memberships.dataobj)).all()

    manifest = json.loads(result.manifest_path.read_text())
    summary = json.loads(result.summary_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["settings"]["source_zero_variance_nan"] is False
    assert summary["loop_entry_gm_flag"] is True
    assert summary["first_thinning_events"] > 0
    assert summary["nonfinite_membership_values"] == 0
    assert summary["source_zero_variance_nan"] is False
    assert summary["zero_lesion_standard_deviation_channels"] == 2
    assert summary["stage_order"][-1] == "propagateThinning"

    # The selected source-faithful default still propagates zero-variance
    # NaNs. The explicit finite override must also reach the *next* GM update
    # so a caller can validate a stabilized iterative branch end to end.
    initial_summary = json.loads(recrop.summary_path.read_text())
    images = np.asarray(nib.load(str(recrop.modalities_path)).dataobj, dtype=np.float32)
    shapes = np.asarray(nib.load(str(recrop.shapes_path)).dataobj, dtype=np.float32)
    transition_inputs = dict(
        centroid_options={
            "mode": "none",
            "modalities": ("T1_MPRAGE", "FLAIR"),
            "lesion_power_exponents": (0.4, 1.2),
            "tissue_power_exponents": np.ones((2, 13), dtype=np.float32),
        },
        channel_weight_options={"outlier_distance": 0.1, "norm_type": "Identity"},
        outlier_energy=0.1,
        lesion_power_exponent=1.2,
    )
    source_faithful_transition = run_lesion_toads_registered_first_iteration(
        captured["initial"],
        images,
        np.asarray(initial_summary["initial_class_centroids"], dtype=np.float32),
        np.asarray(initial_summary["lesion_centroid"], dtype=np.float32),
        np.asarray(initial_summary["black_hole_centroid"], dtype=np.float32),
        shapes,
        np.ones(13, dtype=np.float32),
        np.ones(13, dtype=np.float32),
        **transition_inputs,
    )
    assert np.isnan(source_faithful_transition.membership.membership.objective)
    assert np.isnan(source_faithful_transition.membership.energy_ratio.energy_ratio)
    transition = run_lesion_toads_registered_first_iteration(
        captured["initial"],
        images,
        np.asarray(initial_summary["initial_class_centroids"], dtype=np.float32),
        np.asarray(initial_summary["lesion_centroid"], dtype=np.float32),
        np.asarray(initial_summary["black_hole_centroid"], dtype=np.float32),
        shapes,
        np.ones(13, dtype=np.float32),
        np.ones(13, dtype=np.float32),
        source_zero_variance_nan=False,
        **transition_inputs,
    )
    assert np.isfinite(transition.membership.membership.class_memberships).all()
    assert np.isfinite(transition.membership.membership.objective)
    assert np.isfinite(transition.membership.energy_ratio.energy_ratio)
    shape_atlas = load_lesion_atlas_shape_priors(ATLAS)
    atlas_topology = load_lesion_atlas_topology(ATLAS)
    image_spacing = tuple(
        np.linalg.norm(np.asarray(nib.load(str(recrop.modalities_path)).affine[:3, :3]), axis=0)
    )
    template_spacing = tuple(
        np.linalg.norm(np.asarray(atlas_topology.image.affine[:3, :3]), axis=0)
    )
    branch = run_lesion_toads_registered_first_iteration_branch(
        captured["initial"],
        transition,
        tuple(shape_atlas.shape_maps[item.name] for item in shape_atlas.structures),
        np.asarray(initial_summary["atlas_transform"], dtype=np.float32),
        image_resolutions_mm=image_spacing,
        template_resolutions_mm=template_spacing,
        object_types=STANDARD13_TYPES,
        classification=np.asarray(nib.load(str(recrop.topology_path)).dataobj),
        template_labels=np.asarray([item.label for item in shape_atlas.structures]),
        lesion_class=10,
        iteration=1,
        n_iter_max=3,
        max_distance=0.0,
        algorithm_mode="diffeomorphism",
        relation_limit=3,
        align_shapes=True,
        init_align_iter=1,
    )
    assert branch.registration.registration.transform.shape == (6,)
    assert branch.branch.branch.stop is False
    assert branch.branch.branch.growing is not None
    assert branch.branch.branch.relations is not None

    def unexpected(*args, **kwargs):
        raise AssertionError("valid resume must skip first membership/thinning")

    monkeypatch.setattr(first_stage, "run_lesion_toads_registered_initial_membership", unexpected)
    resumed = first_stage.run_plus_lesion_first_iteration_checkpoint(
        recrop.manifest_path,
        ATLAS,
        tmp_path / "first-iteration",
        iterations=1,
        first_thinning_limit=1000.0,
        source_zero_variance_nan=False,
        resume=True,
    )
    assert resumed.resumed


def test_first_iteration_cli_routes_source_faithful_default(tmp_path: Path, monkeypatch) -> None:
    from pymacruise.cli import main

    captured = {}

    def fake_first(manifest, descriptor, output, **kwargs):
        captured.update(manifest=manifest, descriptor=descriptor, output=output, **kwargs)
        return first_stage.PlusLesionFirstIterationResult(
            tmp_path / "manifest.json", *(tmp_path / f"{name}.nii.gz" for name in (
                "memberships", "segmentation", "skeleton", "availability", "ordering"
            )), tmp_path / "summary.json", True,
        )

    monkeypatch.setattr(first_stage, "run_plus_lesion_first_iteration_checkpoint", fake_first)
    assert main(
        [
            "plus-lesion-first-iteration",
            "--registered-recrop-manifest", "recrop.json",
            "--atlas-descriptor", str(ATLAS),
            "--output-dir", "out",
            "--iterations", "1",
            "--resume",
        ]
    ) == 0
    assert captured["source_zero_variance_nan"] is True
    assert captured["iterations"] == 1
