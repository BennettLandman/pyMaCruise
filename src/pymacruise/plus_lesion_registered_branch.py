"""Restartable first registered post-thinning LesionTOADS branch for plus inputs.

This checkpoint composes the already selected registered initial-membership,
first-transition, atlas-registration, and wrapper branch adapters. It writes a
compact RAS image handoff plus a provenance manifest. It is one synthetic
branch boundary, not a complete LesionTOADS workflow.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from .image_io import load_nifti_ras, nifti_ras_on_reference, save_nifti_ras
from .lesion_atlas_priors import (
    load_lesion_atlas_priors,
    load_lesion_atlas_shape_priors,
    load_lesion_atlas_topology,
)
from .lesion_initial_segmentation_state import initialize_lesion_segmentation_state
from .lesion_registered_first_iteration import run_lesion_toads_registered_first_iteration
from .lesion_registered_first_iteration_branch import run_lesion_toads_registered_first_iteration_branch
from .lesion_registered_initial_membership import (
    LesionInitialMembershipNonfiniteError,
    run_lesion_toads_registered_initial_membership,
)
from .plus_lesion_first_iteration import hydrate_plus_lesion_first_iteration_state
from .lesion_standard13_voxel import STANDARD13_TYPES


@dataclass(frozen=True)
class PlusLesionRegisteredBranchResult:
    manifest_path: Path
    memberships_path: Path
    segmentation_path: Path
    skeleton_path: Path | None
    availability_path: Path | None
    ordering_path: Path | None
    relation_map_path: Path | None
    summary_path: Path
    resumed: bool
    classification_path: Path | None = None
    lesion_subtype_path: Path | None = None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _valid_cached(manifest: dict[str, Any], signature: str, root: Path) -> PlusLesionRegisteredBranchResult | None:
    if manifest.get("status") != "complete" or manifest.get("run_signature") != signature:
        return None
    paths: dict[str, Path] = {}
    for key in (
        "memberships", "membership_diagnostics", "distance_state", "segmentation",
        "skeleton", "availability", "ordering", "relation_map", "summary",
    ):
        record = manifest.get("outputs", {}).get(key)
        if not isinstance(record, dict):
            return None
        path = Path(record.get("path", ""))
        if not path.is_relative_to(root) or not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
            return None
        paths[key] = path
    try:
        images = [nib.load(str(paths[key])) for key in ("memberships", "segmentation", "skeleton", "availability", "ordering", "relation_map")]
        summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return None
    shape = images[0].shape[:3]
    diagnostics = nib.load(str(paths["membership_diagnostics"]))
    distances = nib.load(str(paths["distance_state"]))
    if (
        images[0].shape != (*shape, 17)
        or diagnostics.shape != (*shape, 17)
        or distances.shape != (*shape, 6)
        or images[1].shape != shape
        or images[2].shape != shape
        or images[3].shape != (*shape, 13)
        or images[4].shape != shape
        or images[5].shape != shape
        or any(nib.aff2axcodes(image.affine) != ("R", "A", "S") for image in images)
        or not all(np.isfinite(np.asarray(image.dataobj)).all() for image in images)
        or not np.isfinite(np.asarray(diagnostics.dataobj)).all()
        or not np.isfinite(np.asarray(distances.dataobj)).all()
        or summary.get("spatial_convention") != "RAS"
    ):
        return None
    return PlusLesionRegisteredBranchResult(
        Path(manifest["manifest_path"]), paths["memberships"], paths["segmentation"],
        paths["skeleton"], paths["availability"], paths["ordering"],
        paths["relation_map"], paths["summary"], True,
    )


def run_plus_lesion_registered_branch_checkpoint(
    registered_recrop_manifest_path: str | Path,
    atlas_descriptor_path: str | Path,
    output_dir: str | Path,
    *,
    iterations: int = 99,
    first_thinning_limit: float = 3.0,
    first_thinning_volume_ratio: float = 0.5,
    source_zero_variance_nan: bool = True,
    max_distance: float = 0.0,
    first_iteration_manifest_path: str | Path | None = None,
    resume: bool = False,
) -> PlusLesionRegisteredBranchResult:
    """Run the selected first post-thinning plus update through branch control.

    ``source_zero_variance_nan=True`` preserves the current source-faithful
    default; callers may explicitly choose the finite limiting-similarity
    mode. The selected checkpoint profile disables field correction and uses
    T1/FLAIR centroids with identity channel weighting on this synthetic
    integration seam.
    """
    if isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 1:
        raise ValueError("iterations must be a positive integer")
    if not np.isfinite(first_thinning_limit) or first_thinning_limit < 0:
        raise ValueError("first_thinning_limit must be finite and nonnegative")
    if not np.isfinite(first_thinning_volume_ratio) or not 0 <= first_thinning_volume_ratio <= 1:
        raise ValueError("first_thinning_volume_ratio must be between zero and one")
    if not isinstance(source_zero_variance_nan, bool) or not isinstance(resume, bool):
        raise TypeError("source_zero_variance_nan and resume must be boolean")
    if not np.isfinite(max_distance) or max_distance < 0:
        raise ValueError("max_distance must be finite and nonnegative")
    upstream_path = Path(registered_recrop_manifest_path).expanduser().resolve()
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    for role, path in (("registered-recrop manifest", upstream_path), ("atlas descriptor", descriptor)):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role}: {path}")
    upstream = json.loads(upstream_path.read_text(encoding="utf-8"))
    if upstream.get("status") != "complete" or upstream.get("stage") != "plus-lesion-registered-recrops":
        raise ValueError("input manifest must describe a complete registered recrop stage")
    if upstream.get("inputs", {}).get("atlas_descriptor", {}).get("path") != str(descriptor):
        raise ValueError("registered recrop does not reference this atlas descriptor")
    if upstream["inputs"]["atlas_descriptor"].get("sha256") != _sha256(descriptor):
        raise ValueError("atlas descriptor no longer matches the registered recrop")
    first_manifest_path = (
        None if first_iteration_manifest_path is None
        else Path(first_iteration_manifest_path).expanduser().resolve()
    )
    first_manifest: dict[str, Any] | None = None
    if first_manifest_path is not None:
        if not first_manifest_path.is_file() or first_manifest_path.is_symlink():
            raise FileNotFoundError(f"Missing regular first-iteration checkpoint: {first_manifest_path}")
        first_manifest = json.loads(first_manifest_path.read_text(encoding="utf-8"))
        linked_recrop = first_manifest.get("inputs", {}).get("registered_recrop_manifest", {})
        if (
            first_manifest.get("stage") != "plus-lesion-first-iteration"
            or first_manifest.get("status") != "complete"
            or Path(linked_recrop.get("path", "")).resolve() != upstream_path
            or linked_recrop.get("sha256") != _sha256(upstream_path)
        ):
            raise ValueError("first-iteration checkpoint must be complete and reference this registered recrop")
        first_settings = first_manifest.get("settings", {})
        if (
            float(first_settings.get("first_thinning_limit", -1)) != float(first_thinning_limit)
            or float(first_settings.get("first_thinning_volume_ratio", -1)) != float(first_thinning_volume_ratio)
            or bool(first_settings.get("source_zero_variance_nan")) != source_zero_variance_nan
        ):
            raise ValueError("first-iteration checkpoint settings must match the registered branch request")
        hydrate_plus_lesion_first_iteration_state(first_manifest_path)
    signature_payload = {
        "stage": "plus-lesion-registered-branch",
        "upstream": {"path": str(upstream_path), "sha256": _sha256(upstream_path)},
        "atlas": {"path": str(descriptor), "sha256": _sha256(descriptor)},
        "first_iteration": (
            None if first_manifest_path is None else {
                "path": str(first_manifest_path), "sha256": _sha256(first_manifest_path),
            }
        ),
        "settings": {
            "iterations": iterations,
            "first_thinning_limit": float(first_thinning_limit),
            "first_thinning_volume_ratio": float(first_thinning_volume_ratio),
            "source_zero_variance_nan": source_zero_variance_nan,
            "max_distance": float(max_distance),
        },
    }
    signature = hashlib.sha256(json.dumps(signature_payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "lesion-registered-branch.json"
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Registered-branch checkpoint exists; use resume=True: {manifest_path}")
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old.get("run_signature") != signature:
            raise ValueError("Cannot resume registered branch: inputs or settings differ")
        cached = _valid_cached(old, signature, root)
        if cached is not None:
            return cached
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no registered-branch manifest exists: {manifest_path}")

    records = upstream.get("outputs", {})
    inputs: dict[str, Path] = {}
    for key in ("modalities", "topology", "shapes", "summary"):
        record = records.get(key, {})
        path = Path(record.get("path", "")).resolve()
        if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
            raise ValueError(f"registered recrop {key} output is missing or hash-mismatched")
        inputs[key] = path
    modalities_image = load_nifti_ras(inputs["modalities"])
    topology_image = load_nifti_ras(inputs["topology"])
    shapes_image = load_nifti_ras(inputs["shapes"])
    images = np.asarray(modalities_image.dataobj, dtype=np.float32)
    topology = np.asarray(topology_image.dataobj)
    shapes = np.asarray(shapes_image.dataobj, dtype=np.float32)
    summary = json.loads(inputs["summary"].read_text(encoding="utf-8"))
    if (
        images.ndim != 4 or images.shape[-1] != 2
        or topology.shape != images.shape[:3]
        or shapes.shape != (*images.shape[:3], 13)
        or any(nib.aff2axcodes(image.affine) != ("R", "A", "S") for image in (modalities_image, topology_image, shapes_image))
        or not all(np.allclose(image.affine, modalities_image.affine, rtol=0, atol=1e-4) for image in (topology_image, shapes_image))
        or not np.isfinite(images).all() or not np.isfinite(shapes).all()
        or summary.get("spatial_convention") != "RAS"
    ):
        raise ValueError("registered recrop products must be finite congruent RAS images")
    priors = load_lesion_atlas_priors(descriptor, ("T1", "FLAIR"))
    shape_priors = load_lesion_atlas_shape_priors(descriptor)
    if priors.structures != shape_priors.structures or len(priors.structures) != 13:
        raise ValueError("plus checkpoint requires the selected 13-class atlas order")
    constructor = initialize_lesion_segmentation_state(topology, tuple(item.label for item in priors.structures))
    centers = np.asarray(summary["initial_class_centroids"], dtype=np.float32)
    lesion_centroid = np.asarray(summary["lesion_centroid"], dtype=np.float32)
    black_centroid = np.asarray(summary["black_hole_centroid"], dtype=np.float32)
    if centers.shape != (2, 13) or lesion_centroid.shape != (2,) or black_centroid.shape != (2,):
        raise ValueError("registered recrop lacks valid retained centroid state")

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": "plus-lesion-registered-branch",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "running",
        "spatial_convention": "RAS",
        "settings": signature_payload["settings"],
        "inputs": {
            "registered_recrop_manifest": {"path": str(upstream_path), "sha256": _sha256(upstream_path)},
            "atlas_descriptor": {"path": str(descriptor), "sha256": _sha256(descriptor)},
            **({} if first_manifest_path is None else {
                "first_iteration_manifest": {
                    "path": str(first_manifest_path), "sha256": _sha256(first_manifest_path),
                }
            }),
            **{key: {"path": str(path), "sha256": _sha256(path)} for key, path in inputs.items()},
        },
        "source_order": [
            "computeMemberships:first-GM-enabled",
            "propagateThinning:first",
            "computeCentroids",
            "computeVariances",
            "computeChannelWeights",
            "computeMemberships:post-thinning-GM-enabled",
            "registerShapes",
            "energyRatio",
            "convergence_check",
            "iteration_limit_check",
            "propagateGrowing:if-continuing",
            "computeRelations:if-continuing",
        ],
        "equivalence_note": (
            "Synthetic selected API composition through one registered branch; "
            + ("reuses persisted initial-membership/thinning state; " if first_manifest_path else "recomputes initial membership/thinning; ")
            + "not a full LesionTOADS loop or composite Java match."
        ),
    }
    _json(manifest_path, manifest)
    try:
        if first_manifest_path is None:
            initial = run_lesion_toads_registered_initial_membership(
                images, constructor.segmentation, centers, lesion_centroid, black_centroid,
                shapes, constructor.skeleton, constructor.available, constructor.ordering,
                iterations=iterations,
                first_thinning_limit=first_thinning_limit,
                first_thinning_volume_ratio=first_thinning_volume_ratio,
                source_zero_variance_nan=source_zero_variance_nan,
            )
        else:
            initial = hydrate_plus_lesion_first_iteration_state(first_manifest_path)
        transition = run_lesion_toads_registered_first_iteration(
            initial,
            images,
            centers,
            lesion_centroid,
            black_centroid,
            shapes,
            np.ones(13, dtype=np.float32),
            np.ones(13, dtype=np.float32),
            centroid_options={
                "mode": "none",
                "modalities": ("T1_MPRAGE", "FLAIR"),
                "lesion_power_exponents": (0.4, 1.2),
                "tissue_power_exponents": np.ones((2, 13), dtype=np.float32),
            },
            channel_weight_options={"outlier_distance": 0.1, "norm_type": "Identity"},
            outlier_energy=0.1,
            lesion_power_exponent=1.2,
            source_zero_variance_nan=source_zero_variance_nan,
            atlas_names=tuple(item.name for item in priors.structures),
        )
        if not np.isfinite(transition.membership.membership.objective):
            raise ValueError("post-thinning membership objective is nonfinite; source-faithful zero-variance behavior recorded")
        template_topology = load_lesion_atlas_topology(descriptor)
        image_spacing = tuple(float(v) for v in np.linalg.norm(modalities_image.affine[:3, :3], axis=0))
        template_spacing = tuple(float(v) for v in np.linalg.norm(template_topology.image.affine[:3, :3], axis=0))
        branch_result = run_lesion_toads_registered_first_iteration_branch(
            initial,
            transition,
            tuple(shape_priors.shape_maps[item.name] for item in shape_priors.structures),
            np.asarray(summary["atlas_transform"], dtype=np.float32),
            image_resolutions_mm=image_spacing,
            template_resolutions_mm=template_spacing,
            object_types=STANDARD13_TYPES,
            classification=topology,
            template_labels=np.asarray([item.label for item in priors.structures], dtype=np.int16),
            lesion_class=10,
            iteration=1,
            n_iter_max=max(iterations, 2),
            max_distance=max_distance,
            algorithm_mode="diffeomorphism",
            relation_limit=3,
            align_shapes=True,
            init_align_iter=1,
        )
        branch = branch_result.branch.branch
        if branch.growing is None or branch.relations is None:
            raise RuntimeError("first branch did not produce the expected continuing growth/relation state")
        reference = nib.Nifti1Image(images[..., 0], modalities_image.affine, header=modalities_image.header)
        membership_values = np.concatenate((
            transition.membership.membership.class_memberships,
            transition.membership.membership.lesion_membership[..., None],
            transition.membership.membership.black_hole_membership[..., None],
            transition.membership.membership.outlier_membership[..., None],
            transition.membership.membership.lesion_relation[..., None],
        ), axis=-1)
        membership = transition.membership.membership
        membership_diagnostics = np.concatenate((
            membership.denominator[..., None],
            membership.class_energies,
            membership.lesion_energy[..., None],
            membership.black_hole_energy[..., None],
            membership.running_objective[..., None],
        ), axis=-1).astype(np.float32, copy=False)
        distance_maps = transition.refreshed_distance_maps
        if distance_maps is None:
            raise RuntimeError("GM-enabled post-thinning transition did not refresh distance maps")
        distance_state = np.stack((
            distance_maps.ventricle_mask,
            distance_maps.ventricle,
            distance_maps.gray_matter_mask,
            distance_maps.gray_matter,
            distance_maps.brainstem_mask,
            distance_maps.brainstem,
        ), axis=-1).astype(np.int16, copy=False)
        if membership_diagnostics.shape != (*images.shape[:3], 17):
            raise ValueError("post-thinning membership diagnostics do not have the expected 17 channels")
        if distance_state.shape != (*images.shape[:3], 6):
            raise ValueError("initial distance state does not have the expected 6 channels")
        outputs = {
            "memberships": save_nifti_ras(nifti_ras_on_reference(reference, membership_values, dtype=np.float32), root / "registered_memberships_17ch.nii.gz", dtype=np.float32),
            "membership_diagnostics": save_nifti_ras(nifti_ras_on_reference(reference, membership_diagnostics, dtype=np.float32), root / "registered_membership_diagnostics_17ch.nii.gz", dtype=np.float32),
            "distance_state": save_nifti_ras(nifti_ras_on_reference(reference, distance_state, dtype=np.int16), root / "initial_distance_state_6ch.nii.gz", dtype=np.int16),
            "segmentation": save_nifti_ras(nifti_ras_on_reference(reference, branch.segmentation, dtype=np.int16), root / "registered_segmentation.nii.gz", dtype=np.int16),
            "skeleton": save_nifti_ras(nifti_ras_on_reference(reference, branch.growing.skeleton, dtype=np.int16), root / "registered_skeleton.nii.gz", dtype=np.int16),
            "availability": save_nifti_ras(nifti_ras_on_reference(reference, branch.growing.available.astype(np.uint8), dtype=np.uint8), root / "registered_availability_13ch.nii.gz", dtype=np.uint8),
            "ordering": save_nifti_ras(nifti_ras_on_reference(reference, branch.growing.ordering, dtype=np.float32), root / "registered_ordering.nii.gz", dtype=np.float32),
            "relation_map": save_nifti_ras(nifti_ras_on_reference(reference, branch.relations.relation_map, dtype=np.int16), root / "registered_relation_map.nii.gz", dtype=np.int16),
        }
        summary_path = root / "registered_branch_summary.json"
        summary_output = {
            "spatial_convention": "RAS",
            "class_names": [item.name for item in priors.structures],
            "transition_objective": float(transition.membership.membership.objective),
            "energy_ratio": float(transition.membership.energy_ratio.energy_ratio),
            "energy_ratio_state": {
                key: float(getattr(transition.membership.energy_ratio, key))
                for key in (
                    "previous_energy", "current_energy", "difference", "numerator",
                    "denominator", "quotient", "energy_ratio",
                )
            },
            "refreshed_centroids": {
                "class_centroids": transition.statistics.statistics.centroids.class_centroids.tolist(),
                "lesion_centroids": transition.statistics.statistics.centroids.lesion_centroids.tolist(),
                "black_hole_centroids": transition.statistics.statistics.centroids.black_hole_centroids.tolist(),
            },
            "registration_transform": branch_result.registration.registration.transform.tolist(),
            "registration_rotation": branch_result.registration.registration.rotation.tolist(),
            "registration_shape_bounds": [
                {
                    "lower_xyz": list(bounds.lower_xyz),
                    "upper_xyz": list(bounds.upper_xyz),
                    "support_count": bounds.support_count,
                }
                for bounds in branch_result.registration.registration.shape_bounds
            ],
            "stop": branch.stop,
            "stop_reason": branch.stop_reason,
            "growing_events": len(branch.growing.events),
            "relation_count": branch.relations.relation_count,
            "relation_diffusion_iterations": branch.relations.diffusion_iterations,
            "relation_membership": branch.relations.relation_membership.astype(np.uint8).tolist(),
            "class_volume": branch.growing.class_volume.tolist(),
            "stage_order": list(branch_result.stage_order),
            "source_zero_variance_nan": source_zero_variance_nan,
            "equivalence_note": manifest["equivalence_note"],
        }
        _json(summary_path, summary_output)
        outputs["summary"] = summary_path
        manifest["outputs"] = {key: {"path": str(path), "sha256": _sha256(path)} for key, path in outputs.items()}
        manifest["status"] = "complete"
        _json(manifest_path, manifest)
    except LesionInitialMembershipNonfiniteError as error:
        manifest["status"] = "failed"
        manifest["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
            "source_zero_variance_nan": True,
            "initial_objective_finite": False,
            "nonfinite_initial_membership_values": int(np.count_nonzero(~np.isfinite(error.membership.class_memberships))),
            "nonfinite_initial_energy_values": int(
                np.count_nonzero(~np.isfinite(error.membership.class_energies))
                + np.count_nonzero(~np.isfinite(error.membership.lesion_energy))
                + np.count_nonzero(~np.isfinite(error.membership.black_hole_energy))
                + np.count_nonzero(~np.isfinite(error.membership.running_objective))
                + np.count_nonzero(~np.isfinite(error.membership.denominator))
            ),
            "zero_class_standard_deviation_channels": int(
                np.count_nonzero(error.variances.class_standard_deviation == np.float32(0.0))
            ),
            "zero_lesion_standard_deviation_channels": int(
                np.count_nonzero(error.variances.lesion_standard_deviation == np.float32(0.0))
            ),
        }
        _json(manifest_path, manifest)
        raise
    except Exception as error:
        manifest["status"] = "failed"
        manifest["failure"] = {"type": type(error).__name__, "message": str(error)}
        _json(manifest_path, manifest)
        raise
    return PlusLesionRegisteredBranchResult(
        manifest_path, outputs["memberships"], outputs["segmentation"], outputs["skeleton"],
        outputs["availability"], outputs["ordering"], outputs["relation_map"], summary_path, False,
    )
