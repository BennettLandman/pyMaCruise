"""Restartable first registered LesionTOADS membership/thinning checkpoint."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import nibabel as nib
import numpy as np

from .image_io import load_nifti_ras, nifti_ras_on_reference, save_nifti_ras
from .lesion_atlas_priors import load_lesion_atlas_priors
from .lesion_initial_segmentation_state import initialize_lesion_segmentation_state
from .lesion_initial_loop_entry import seed_lesion_toads_initial_loop
from .lesion_initial_thinning_composition import LesionInitialThinningComposition
from .lesion_growth import (
    LesionToadsThinningInitialization,
    LesionToadsThinningPass,
    LesionToadsThinningQueueRun,
)
from .lesion_registered_initial_membership import (
    LesionInitialDistanceMaps,
    LesionRegisteredInitialMembership,
)
from .json_numbers import json_number
from .lesion_relations import LesionRelationList
from .lesion_standard13_volume import Standard13VolumeResult
from .lesion_registered_initial_membership import (
    LesionInitialMembershipNonfiniteError,
    run_lesion_toads_registered_initial_membership,
)


@dataclass(frozen=True)
class PlusLesionFirstIterationResult:
    """First membership outputs and thinned state checkpoint."""

    manifest_path: Path
    memberships_path: Path
    segmentation_path: Path
    skeleton_path: Path
    availability_path: Path
    ordering_path: Path
    summary_path: Path
    resumed: bool


def hydrate_plus_lesion_first_iteration_state(
    first_iteration_manifest_path: str | Path,
) -> LesionRegisteredInitialMembership:
    """Reconstruct P41i membership/thinning state from its hashed outputs."""
    manifest_path = Path(first_iteration_manifest_path).expanduser().resolve()
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise FileNotFoundError(f"Missing regular first-iteration manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete" or manifest.get("stage") != "plus-lesion-first-iteration":
        raise ValueError("input must be a complete plus first-iteration checkpoint")
    outputs: dict[str, Path] = {}
    for name in (
        "memberships", "membership_diagnostics", "distance_state", "relation_map",
        "segmentation", "skeleton", "availability", "ordering", "summary",
    ):
        record = manifest.get("outputs", {}).get(name, {})
        path = Path(record.get("path", "")).resolve()
        if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
            raise ValueError(f"first-iteration {name} output is missing or hash-mismatched")
        outputs[name] = path
    summary = json.loads(outputs["summary"].read_text(encoding="utf-8"))
    if summary.get("spatial_convention") != "RAS":
        raise ValueError("first-iteration summary must declare RAS")
    images = {name: load_nifti_ras(path) for name, path in outputs.items() if name != "summary"}
    reference = images["memberships"]
    spatial = reference.shape[:3]
    for name, image in images.items():
        if (
            nib.aff2axcodes(image.affine) != ("R", "A", "S")
            or image.shape[:3] != spatial
            or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4)
        ):
            raise ValueError(f"first-iteration {name} output is not on the shared RAS grid")
    packed = np.asarray(images["memberships"].dataobj, dtype=np.float32)
    diagnostics = np.asarray(images["membership_diagnostics"].dataobj, dtype=np.float32)
    distances = np.asarray(images["distance_state"].dataobj, dtype=np.int16)
    if packed.shape != (*spatial, 17) or diagnostics.shape != (*spatial, 17) or distances.shape != (*spatial, 6):
        raise ValueError("first-iteration checkpoint state has invalid channel counts")
    diagnostic_nonfinite = int(np.count_nonzero(~np.isfinite(diagnostics)))
    if (
        not np.isfinite(packed).all()
        or diagnostic_nonfinite != int(summary.get("nonfinite_membership_energy_values", -1))
        or (diagnostic_nonfinite and not manifest.get("settings", {}).get("source_zero_variance_nan"))
    ):
        raise ValueError("first-iteration membership state contains unreported nonfinite values")
    membership = Standard13VolumeResult(
        class_memberships=packed[..., :13],
        lesion_membership=packed[..., 13],
        black_hole_membership=packed[..., 14],
        outlier_membership=packed[..., 15],
        lesion_relation=packed[..., 16],
        denominator=diagnostics[..., 0],
        class_energies=diagnostics[..., 1:14],
        lesion_energy=diagnostics[..., 14],
        black_hole_energy=diagnostics[..., 15],
        objective=np.float32(float(summary["initial_energy"])),
        running_objective=diagnostics[..., 16],
    )
    loop_entry = seed_lesion_toads_initial_loop(
        membership, iterations=int(manifest["settings"]["iterations"])
    )
    thinning_run = LesionToadsThinningQueueRun(
        segmentation=np.asarray(images["segmentation"].dataobj, dtype=np.int16),
        skeleton=np.asarray(images["skeleton"].dataobj, dtype=np.int16),
        available=np.asarray(images["availability"].dataobj, dtype=np.bool_),
        ordering=np.asarray(images["ordering"].dataobj, dtype=np.float32),
        class_volume=np.asarray(summary["first_thinning_class_volume"], dtype=np.float32),
        events=(),
        queue=None,
    )
    thinning_initialization = LesionToadsThinningInitialization(
        skeleton=thinning_run.skeleton,
        available=thinning_run.available,
        class_volume=thinning_run.class_volume,
        boundary_counts=np.zeros(spatial, dtype=np.int32),
        inserted_items=(),
        queue=None,
    )
    thinning = LesionToadsThinningPass(
        initialization=thinning_initialization,
        minimum_class_volume=np.zeros(13, dtype=np.float32),
        run=thinning_run,
    )
    first_thinning = LesionInitialThinningComposition(
        loop_entry=loop_entry,
        thinning=thinning,
        stage_order=tuple(summary["stage_order"]),
    )
    relations = LesionRelationList(
        np.asarray(images["relation_map"].dataobj, dtype=np.int16),
        np.asarray(summary["relation_membership"], dtype=np.bool_),
        diffusion_iterations=int(summary["relation_diffusion_iterations"]),
    )
    distance_maps = LesionInitialDistanceMaps(
        distances[..., 0].astype(np.bool_), distances[..., 1].astype(np.int8),
        distances[..., 2].astype(np.bool_), distances[..., 3].astype(np.int8),
        distances[..., 4].astype(np.bool_), distances[..., 5].astype(np.int8),
    )
    return LesionRegisteredInitialMembership(
        relations=relations,
        variances=SimpleNamespace(class_standard_deviation=None, lesion_standard_deviation=None),
        channel_weights=np.empty((0,), dtype=np.float32),
        distance_maps=distance_maps,
        initial_class_memberships=membership.class_memberships,
        initial_lesion_membership=membership.lesion_membership,
        initial_black_hole_membership=membership.black_hole_membership,
        initial_outlier_membership=membership.outlier_membership,
        initial_lesion_relation=membership.lesion_relation,
        membership=membership,
        loop_entry=loop_entry,
        first_thinning=first_thinning,
        stage_order=tuple(summary["stage_order"]),
        excluded_wrapper_stages=(
            "atlas registration supplied by registered recrop",
            "inhomogeneity correction disabled in selected plus profile",
        ),
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _signature(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _cached(manifest: dict[str, Any], signature: str, root: Path) -> PlusLesionFirstIterationResult | None:
    if manifest.get("status") != "complete" or manifest.get("run_signature") != signature:
        return None
    paths: dict[str, Path] = {}
    records = manifest.get("outputs", {})
    for name in (
        "memberships", "membership_diagnostics", "distance_state", "relation_map",
        "segmentation", "skeleton", "availability", "ordering", "summary",
    ):
        record = records.get(name)
        if not isinstance(record, dict):
            return None
        path = Path(record.get("path", ""))
        if not path.is_relative_to(root) or not path.is_file() or path.is_symlink():
            return None
        if _sha256(path) != record.get("sha256"):
            return None
        paths[name] = path
    try:
        memberships = nib.load(str(paths["memberships"]))
        diagnostics = nib.load(str(paths["membership_diagnostics"]))
        distances = nib.load(str(paths["distance_state"]))
        relations = nib.load(str(paths["relation_map"]))
        segmentation = nib.load(str(paths["segmentation"]))
        skeleton = nib.load(str(paths["skeleton"]))
        availability = nib.load(str(paths["availability"]))
        ordering = nib.load(str(paths["ordering"]))
        summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return None
    spatial = memberships.shape[:3]
    diagnostics_array = np.asarray(diagnostics.dataobj)
    diagnostics_nonfinite = int(np.count_nonzero(~np.isfinite(diagnostics_array)))
    if (
        len(memberships.shape) != 4
        or memberships.shape[-1] != 17
        or diagnostics.shape != (*spatial, 17)
        or distances.shape != (*spatial, 6)
        or relations.shape != spatial
        or segmentation.shape != spatial
        or skeleton.shape != spatial
        or ordering.shape != spatial
        or len(availability.shape) != 4
        or availability.shape != (*spatial, 13)
        or any(nib.aff2axcodes(image.affine) != ("R", "A", "S") for image in (memberships, diagnostics, distances, relations, segmentation, skeleton, availability, ordering))
        or not np.isfinite(np.asarray(memberships.dataobj)).all()
        or diagnostics_nonfinite != int(summary.get("nonfinite_membership_energy_values", -1))
        or (diagnostics_nonfinite and not manifest.get("settings", {}).get("source_zero_variance_nan"))
        or summary.get("spatial_convention") != "RAS"
    ):
        return None
    return PlusLesionFirstIterationResult(
        Path(manifest["manifest_path"]), paths["memberships"], paths["segmentation"],
        paths["skeleton"], paths["availability"], paths["ordering"], paths["summary"], True,
    )


def run_plus_lesion_first_iteration_checkpoint(
    registered_recrop_manifest_path: str | Path,
    atlas_descriptor_path: str | Path,
    output_dir: str | Path,
    *,
    iterations: int = 99,
    first_thinning_limit: float = 3.0,
    first_thinning_volume_ratio: float = 0.5,
    source_zero_variance_nan: bool = True,
    resume: bool = False,
) -> PlusLesionFirstIterationResult:
    """Compute registered initial memberships and the first thinning handoff.

    The preceding P41h manifest supplies the cropped T1/FLAIR image, atlas
    labels/shapes, and retained pre-registration class centroids. The selected
    standard-13 adapter computes the initial membership volume with GMflag
    enabled, then creates the source loop-entry state and first thinning pass.
    """
    if isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 1:
        raise ValueError("iterations must be a positive integer")
    if not np.isfinite(first_thinning_limit) or first_thinning_limit < 0:
        raise ValueError("first_thinning_limit must be finite and nonnegative")
    if not np.isfinite(first_thinning_volume_ratio) or not 0 <= first_thinning_volume_ratio <= 1:
        raise ValueError("first_thinning_volume_ratio must be between zero and one")
    if not isinstance(source_zero_variance_nan, bool) or not isinstance(resume, bool):
        raise TypeError("source_zero_variance_nan and resume must be boolean")
    upstream_manifest_path = Path(registered_recrop_manifest_path).expanduser().resolve()
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    for role, path in (("registered-recrop manifest", upstream_manifest_path), ("atlas descriptor", descriptor)):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role}: {path}")
    upstream = json.loads(upstream_manifest_path.read_text(encoding="utf-8"))
    if upstream.get("status") != "complete" or upstream.get("stage") != "plus-lesion-registered-recrops":
        raise ValueError("input manifest must describe a complete registered recrop stage")
    inputs = upstream.get("inputs", {})
    descriptor_input = inputs.get("atlas_descriptor", {})
    if descriptor_input.get("path") != str(descriptor) or descriptor_input.get("sha256") != _sha256(descriptor):
        raise ValueError("registered recrop does not reference this atlas descriptor")
    signature = _signature({
        "stage": "plus-lesion-first-iteration",
        "registered_recrop_manifest": {"path": str(upstream_manifest_path), "sha256": _sha256(upstream_manifest_path)},
        "atlas_descriptor": {"path": str(descriptor), "sha256": _sha256(descriptor)},
        "settings": {
            "iterations": iterations,
            "first_thinning_limit": float(first_thinning_limit),
            "first_thinning_volume_ratio": float(first_thinning_volume_ratio),
            "source_zero_variance_nan": source_zero_variance_nan,
        },
    })
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "lesion-first-iteration.json"
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"First-iteration checkpoint exists; use resume=True: {manifest_path}")
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old.get("run_signature") != signature:
            raise ValueError("Cannot resume first iteration: inputs or settings differ")
        cached = _cached(old, signature, root)
        if cached is not None:
            return cached
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no first-iteration manifest exists: {manifest_path}")

    output_records = upstream.get("outputs", {})
    paths: dict[str, Path] = {}
    for key in ("modalities", "topology", "shapes", "summary"):
        record = output_records.get(key, {})
        path = Path(record.get("path", "")).resolve()
        if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
            raise ValueError(f"registered recrop {key} output is missing or hash-mismatched")
        paths[key] = path
    modalities_image = load_nifti_ras(paths["modalities"])
    topology_image = load_nifti_ras(paths["topology"])
    shapes_image = load_nifti_ras(paths["shapes"])
    summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    images = np.asarray(modalities_image.dataobj, dtype=np.float32)
    topology = np.asarray(topology_image.dataobj)
    shapes = np.asarray(shapes_image.dataobj, dtype=np.float32)
    if (
        images.ndim != 4 or images.shape[-1] != 2
        or topology.shape != images.shape[:3]
        or shapes.shape != (*images.shape[:3], 13)
        or any(nib.aff2axcodes(image.affine) != ("R", "A", "S") for image in (modalities_image, topology_image, shapes_image))
        or any(not np.allclose(image.affine, modalities_image.affine, rtol=0, atol=1e-4) for image in (topology_image, shapes_image))
        or not np.isfinite(images).all() or not np.isfinite(shapes).all()
        or summary.get("spatial_convention") != "RAS"
    ):
        raise ValueError("registered recrop outputs must be finite, congruent RAS images")
    from .lesion_standard13_voxel import STANDARD13_NAMES

    priors = load_lesion_atlas_priors(descriptor, ("T1", "FLAIR"))
    if tuple(item.name for item in priors.structures) != tuple(STANDARD13_NAMES):
        raise ValueError("atlas class order must match the selected standard-13 membership adapter")
    state = initialize_lesion_segmentation_state(
        topology, tuple(item.label for item in priors.structures)
    )
    centers = np.asarray(summary.get("initial_class_centroids"), dtype=np.float32)
    lesion_centroid = np.asarray(summary.get("lesion_centroid"), dtype=np.float32)
    black_hole_centroid = np.asarray(summary.get("black_hole_centroid"), dtype=np.float32)
    if centers.shape != (2, 13) or lesion_centroid.shape != (2,) or black_hole_centroid.shape != (2,):
        raise ValueError("registered recrop summary lacks expected centroid handoff arrays")

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": "plus-lesion-first-iteration",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "running",
        "spatial_convention": "RAS",
        "source_order": [
            "LesionToadSegmentation.constructor:registeredRecropGrid",
            "setGMflag",
            "computeMemberships:initial",
            "Eprev=initialEnergy;Eratio=1",
            "initThinningLabels",
            "propagateThinning",
        ],
        "settings": {
            "iterations": iterations,
            "first_thinning_limit": float(first_thinning_limit),
            "first_thinning_volume_ratio": float(first_thinning_volume_ratio),
            "source_zero_variance_nan": source_zero_variance_nan,
        },
        "inputs": {
            "registered_recrop_manifest": {"path": str(upstream_manifest_path), "sha256": _sha256(upstream_manifest_path)},
            "atlas_descriptor": {"path": str(descriptor), "sha256": _sha256(descriptor)},
        },
        "equivalence_note": "Source-ordered selected synthetic initial membership and first thinning composition; later LesionTOADS iterations and whole-workflow equivalence remain open.",
    }
    _write_json(manifest_path, manifest)
    try:
        result = run_lesion_toads_registered_initial_membership(
            images,
            state.segmentation,
            centers,
            lesion_centroid,
            black_hole_centroid,
            shapes,
            state.skeleton,
            state.available,
            state.ordering,
            iterations=iterations,
            first_thinning_limit=float(first_thinning_limit),
            first_thinning_volume_ratio=float(first_thinning_volume_ratio),
            source_zero_variance_nan=source_zero_variance_nan,
        )
    except LesionInitialMembershipNonfiniteError as error:
        partial_memberships = np.concatenate(
            (
                error.membership.class_memberships,
                error.membership.lesion_membership[..., None],
                error.membership.black_hole_membership[..., None],
                error.membership.outlier_membership[..., None],
                error.membership.lesion_relation[..., None],
            ),
            axis=-1,
        )
        manifest["status"] = "failed"
        manifest["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
            "source_zero_variance_nan": source_zero_variance_nan,
            "initial_objective_finite": False,
            "nonfinite_initial_membership_values": int(
                np.count_nonzero(~np.isfinite(partial_memberships))
            ),
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
        _write_json(manifest_path, manifest)
        raise
    except ValueError as error:
        manifest["status"] = "failed"
        manifest["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
            "source_zero_variance_nan": source_zero_variance_nan,
        }
        _write_json(manifest_path, manifest)
        raise
    reference = nib.Nifti1Image(
        images[..., 0], modalities_image.affine, header=modalities_image.header
    )
    memberships = np.concatenate(
        (
            result.membership.class_memberships,
            result.membership.lesion_membership[..., None],
            result.membership.black_hole_membership[..., None],
            result.membership.outlier_membership[..., None],
            result.membership.lesion_relation[..., None],
        ),
        axis=-1,
    )
    memberships_path = save_nifti_ras(
        nifti_ras_on_reference(reference, memberships, dtype=np.float32),
        root / "initial_memberships_17ch.nii.gz", dtype=np.float32,
    )
    membership = result.membership
    diagnostics = np.concatenate((
        membership.denominator[..., None], membership.class_energies,
        membership.lesion_energy[..., None], membership.black_hole_energy[..., None],
        membership.running_objective[..., None],
    ), axis=-1).astype(np.float32, copy=False)
    diagnostics_path = save_nifti_ras(
        nifti_ras_on_reference(reference, diagnostics, dtype=np.float32),
        root / "initial_membership_diagnostics_17ch.nii.gz", dtype=np.float32,
    )
    distance_maps = result.distance_maps
    distance_state = np.stack((
        distance_maps.ventricle_mask, distance_maps.ventricle,
        distance_maps.gray_matter_mask, distance_maps.gray_matter,
        distance_maps.brainstem_mask, distance_maps.brainstem,
    ), axis=-1).astype(np.int16, copy=False)
    distance_state_path = save_nifti_ras(
        nifti_ras_on_reference(reference, distance_state, dtype=np.int16),
        root / "initial_distance_state_6ch.nii.gz", dtype=np.int16,
    )
    relation_map_path = save_nifti_ras(
        nifti_ras_on_reference(reference, result.relations.relation_map, dtype=np.int16),
        root / "initial_relation_map.nii.gz", dtype=np.int16,
    )
    thinning = result.first_thinning.thinning.run
    segmentation_path = save_nifti_ras(
        nifti_ras_on_reference(reference, thinning.segmentation, dtype=np.int16),
        root / "first_thinning_segmentation.nii.gz", dtype=np.int16,
    )
    skeleton_path = save_nifti_ras(
        nifti_ras_on_reference(reference, thinning.skeleton, dtype=np.int16),
        root / "first_thinning_skeleton.nii.gz", dtype=np.int16,
    )
    availability_path = save_nifti_ras(
        nifti_ras_on_reference(reference, thinning.available.astype(np.uint8), dtype=np.uint8),
        root / "first_thinning_available_13ch.nii.gz", dtype=np.uint8,
    )
    ordering_path = save_nifti_ras(
        nifti_ras_on_reference(reference, thinning.ordering, dtype=np.float32),
        root / "first_thinning_ordering.nii.gz", dtype=np.float32,
    )
    summary_path = root / "first_iteration_summary.json"
    nan_count = int(np.count_nonzero(~np.isfinite(memberships)))
    summary_output = {
        "spatial_convention": "RAS",
        "class_names": list(priors.structures[i].name for i in range(13)),
        "membership_channels": [*list(priors.structures[i].name for i in range(13)), "lesion", "black_hole", "outlier", "lesion_relation"],
        "initial_energy": json_number(result.membership.objective),
        "loop_entry_energy_ratio": json_number(result.loop_entry.energy_ratio),
        "loop_entry_previous_energy": json_number(result.loop_entry.previous_energy),
        "loop_entry_gm_flag": bool(result.loop_entry.gm_flag),
        "first_thinning_events": len(thinning.events),
        "first_thinning_class_volume": thinning.class_volume.tolist(),
        "relation_count": result.relations.relation_count,
        "relation_diffusion_iterations": result.relations.diffusion_iterations,
        "relation_membership": result.relations.relation_membership.astype(np.uint8).tolist(),
        "nonfinite_membership_values": nan_count,
        "nonfinite_membership_energy_values": int(np.count_nonzero(~np.isfinite(diagnostics))),
        "source_zero_variance_nan": source_zero_variance_nan,
        "zero_class_standard_deviation_channels": int(
            np.count_nonzero(result.variances.class_standard_deviation == np.float32(0.0))
        ),
        "zero_lesion_standard_deviation_channels": int(
            np.count_nonzero(result.variances.lesion_standard_deviation == np.float32(0.0))
        ),
        "zero_class_standard_deviation_channels": int(
            np.count_nonzero(result.variances.class_standard_deviation == np.float32(0.0))
        ),
        "zero_lesion_standard_deviation_channels": int(
            np.count_nonzero(result.variances.lesion_standard_deviation == np.float32(0.0))
        ),
        "stage_order": list(result.stage_order),
        "equivalence_note": manifest["equivalence_note"],
    }
    _write_json(summary_path, summary_output)
    manifest["outputs"] = {
        name: {"path": str(path), "sha256": _sha256(path)}
        for name, path in (
            ("memberships", memberships_path), ("membership_diagnostics", diagnostics_path),
            ("distance_state", distance_state_path), ("relation_map", relation_map_path),
            ("segmentation", segmentation_path),
            ("skeleton", skeleton_path), ("availability", availability_path),
            ("ordering", ordering_path), ("summary", summary_path),
        )
    }
    manifest["status"] = "complete"
    _write_json(manifest_path, manifest)
    return PlusLesionFirstIterationResult(
        manifest_path, memberships_path, segmentation_path, skeleton_path,
        availability_path, ordering_path, summary_path, False,
    )
