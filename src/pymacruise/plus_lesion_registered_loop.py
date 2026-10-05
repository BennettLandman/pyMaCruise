"""Continue a persisted synthetic plus branch with the selected P39 loop API.

This adapter reconstructs the in-memory branch state written by P41k. It
does not yet write the next loop boundary; the output remains an API result so
the composition can be compared directly with the existing P39 loop tests.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np

from .image_io import load_nifti_ras, nifti_ras_on_reference, save_nifti_ras
from .lesion_atlas_priors import (
    load_lesion_atlas_priors,
    load_lesion_atlas_shape_priors,
    load_lesion_atlas_topology,
)
from .lesion_descriptor_loop import (
    LesionToadsDescriptorLoopOptions,
    LesionToadsDescriptorLoopResult,
    run_lesion_toads_descriptor_loop,
)
from .lesion_energy_ratio import LesionEnergyRatio
from .lesion_growth import LesionToadsGrowingPass
from .lesion_next_iteration_membership import LesionNextIterationMembership
from .lesion_post_registration import LesionToadsPostRegistrationBranch
from .lesion_registration_handoff import LesionRegistrationHandoff
from .lesion_relations import LesionRelationList
from .lesion_standard13_volume import Standard13VolumeResult
from .toads_shapes import ToadShapeBounds
from .lesion_standard13_voxel import STANDARD13_TYPES
from .lesion_registered_initial_membership import LesionInitialDistanceMaps


class RegisteredBranchStateError(ValueError):
    """Raised when a plus checkpoint cannot be safely hydrated."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _load_complete_checkpoint(manifest_path: str | Path) -> tuple[dict, dict[str, Path]]:
    manifest_file = Path(manifest_path).expanduser().resolve()
    if not manifest_file.is_file() or manifest_file.is_symlink():
        raise FileNotFoundError(f"Missing regular registered-branch manifest: {manifest_file}")
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete" or manifest.get("stage") not in {
        "plus-lesion-registered-branch", "plus-lesion-registered-loop",
    }:
        raise RegisteredBranchStateError("input must be a complete plus registered-branch or loop checkpoint")
    outputs: dict[str, Path] = {}
    for name, record in manifest.get("outputs", {}).items():
        if not isinstance(record, dict):
            continue
        path = Path(record.get("path", "")).resolve()
        if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
            raise RegisteredBranchStateError(f"registered-branch output {name!r} is missing or hash-mismatched")
        outputs[name] = path
    required = {
        "memberships", "membership_diagnostics", "distance_state", "segmentation",
        "skeleton", "availability", "ordering", "relation_map", "summary",
    }
    missing = required - outputs.keys()
    if missing:
        raise RegisteredBranchStateError(f"registered-branch checkpoint lacks outputs: {sorted(missing)}")
    return manifest, outputs


def hydrate_plus_lesion_registered_branch(
    manifest_path: str | Path,
) -> tuple[LesionToadsPostRegistrationBranch, LesionRegistrationHandoff, np.ndarray]:
    """Reconstruct P41k branch objects and six distance maps from hashed files."""
    manifest, outputs = _load_complete_checkpoint(manifest_path)
    summary = json.loads(outputs["summary"].read_text(encoding="utf-8"))
    if summary.get("spatial_convention") != "RAS":
        raise RegisteredBranchStateError("registered-branch summary must declare RAS")
    images = {
        name: load_nifti_ras(path)
        for name, path in outputs.items()
        if name not in {"summary"}
    }
    reference = images["memberships"]
    spatial = reference.shape[:3]
    for name, image in images.items():
        if image.shape[:3] != spatial or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4):
            raise RegisteredBranchStateError(f"registered-branch output {name!r} is not on the shared RAS grid")
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise RegisteredBranchStateError(f"registered-branch output {name!r} is not RAS")

    packed = np.asarray(images["memberships"].dataobj, dtype=np.float32)
    diagnostics = np.asarray(images["membership_diagnostics"].dataobj, dtype=np.float32)
    distance_state = np.asarray(images["distance_state"].dataobj, dtype=np.int16)
    if packed.shape != (*spatial, 17) or diagnostics.shape != (*spatial, 17) or distance_state.shape != (*spatial, 6):
        raise RegisteredBranchStateError("registered-branch state channel shapes are invalid")
    if not all(np.isfinite(array).all() for array in (packed, diagnostics)):
        raise RegisteredBranchStateError("registered-branch membership state contains nonfinite values")

    ratio_state = summary["energy_ratio_state"]
    ratio = LesionEnergyRatio(**{key: np.float32(value) for key, value in ratio_state.items()})
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
        objective=np.float32(summary["transition_objective"]),
        running_objective=diagnostics[..., 16],
    )
    membership_update = LesionNextIterationMembership(membership, ratio, tuple(summary["stage_order"]))

    bounds = tuple(
        ToadShapeBounds(
            tuple(int(value) for value in record["lower_xyz"]),
            tuple(int(value) for value in record["upper_xyz"]),
            int(record["support_count"]),
            None,
        )
        for record in summary["registration_shape_bounds"]
    )
    registration = LesionRegistrationHandoff(
        memberships=membership.class_memberships,
        optimization=None,
        transform=np.asarray(summary["registration_transform"], dtype=np.float32),
        rotation=np.asarray(summary["registration_rotation"], dtype=np.float32),
        shape_bounds=bounds,
    )
    growing = LesionToadsGrowingPass(
        segmentation=np.asarray(images["segmentation"].dataobj, dtype=np.int16),
        skeleton=np.asarray(images["skeleton"].dataobj, dtype=np.int16),
        available=np.asarray(images["availability"].dataobj, dtype=np.bool_),
        ordering=np.asarray(images["ordering"].dataobj, dtype=np.float32),
        class_volume=np.asarray(summary["class_volume"], dtype=np.float32),
        events=(),
        remaining_queue=(),
    )
    relations = LesionRelationList(
        np.asarray(images["relation_map"].dataobj, dtype=np.int16),
        np.asarray(summary["relation_membership"], dtype=np.bool_),
        diffusion_iterations=int(summary.get("relation_diffusion_iterations", 0)),
    )
    branch = LesionToadsPostRegistrationBranch(
        membership_result=membership_update,
        registration_handoff=registration,
        energy_ratio=np.float32(ratio.energy_ratio),
        convergence_stop=False,
        iteration_limit_stop=False,
        stop=False,
        stop_reason=None,
        run_growing=True,
        growing=growing,
        segmentation=growing.segmentation,
        final_classification=None,
        relations=relations,
        stage_order=tuple(summary["stage_order"]),
    )
    return branch, registration, distance_state


def continue_plus_lesion_registered_loop(
    registered_branch_manifest: str | Path,
    atlas_descriptor_path: str | Path,
    *,
    n_iter_max: int = 99,
    thinning_limit: float = 3.0,
    source_zero_variance_nan: bool | None = None,
    max_passes: int | None = None,
) -> LesionToadsDescriptorLoopResult:
    """Hydrate P41k and run the next selected registered P39 iteration.

    The output is an in-memory selected-loop result. Persisting the next loop
    boundary and adding a public CLI are subsequent integration work.
    """
    if isinstance(n_iter_max, bool) or not isinstance(n_iter_max, int) or n_iter_max < 2:
        raise ValueError("n_iter_max must be an integer of at least 2")
    manifest, outputs = _load_complete_checkpoint(registered_branch_manifest)
    if source_zero_variance_nan is None:
        source_zero_variance_nan = bool(manifest["settings"]["source_zero_variance_nan"])
    if not isinstance(source_zero_variance_nan, bool):
        raise TypeError("source_zero_variance_nan must be boolean")
    branch, registration, distance_state = hydrate_plus_lesion_registered_branch(registered_branch_manifest)
    upstream = manifest["inputs"]["registered_recrop_manifest"]
    recrop_manifest_path = Path(upstream["path"]).resolve()
    if not recrop_manifest_path.is_file() or _sha256(recrop_manifest_path) != upstream["sha256"]:
        raise RegisteredBranchStateError("registered-recrop manifest is missing or hash-mismatched")
    recrop = json.loads(recrop_manifest_path.read_text(encoding="utf-8"))
    recorded_atlas = manifest["inputs"]["atlas_descriptor"]
    atlas_record_path = Path(recorded_atlas["path"]).resolve()
    atlas_path = Path(atlas_descriptor_path).expanduser().resolve()
    if atlas_path != atlas_record_path or _sha256(atlas_path) != recorded_atlas["sha256"]:
        raise RegisteredBranchStateError("atlas descriptor differs from the registered-branch input")
    input_paths: dict[str, Path] = {}
    for name in ("modalities", "topology"):
        record = recrop["outputs"][name]
        path = Path(record["path"]).resolve()
        if not path.is_file() or path.is_symlink() or _sha256(path) != record["sha256"]:
            raise RegisteredBranchStateError(f"registered-recrop {name} image is missing or hash-mismatched")
        input_paths[name] = path
    image_path = input_paths["modalities"]
    topology_path = input_paths["topology"]
    image = load_nifti_ras(image_path)
    topology = load_nifti_ras(topology_path)
    if (
        len(image.shape) != 4 or image.shape[-1] != 2 or topology.shape != image.shape[:3]
        or nib.aff2axcodes(image.affine) != ("R", "A", "S")
        or nib.aff2axcodes(topology.affine) != ("R", "A", "S")
        or not np.allclose(image.affine, topology.affine, rtol=0, atol=1e-4)
    ):
        raise RegisteredBranchStateError("registered-recrop modalities and topology must share a RAS grid")
    images = np.asarray(image.dataobj, dtype=np.float32)
    priors = load_lesion_atlas_priors(atlas_path, ("T1", "FLAIR"))
    atlas_shapes = load_lesion_atlas_shape_priors(atlas_path)
    atlas_topology = load_lesion_atlas_topology(atlas_path)
    if len(priors.structures) != 13 or priors.structures != atlas_shapes.structures:
        raise RegisteredBranchStateError("loop continuation requires the selected 13-class atlas order")
    summary = json.loads(outputs["summary"].read_text(encoding="utf-8"))
    starting_iteration = int(summary.get("completed_iteration", 1))
    if n_iter_max <= starting_iteration:
        raise ValueError("n_iter_max must exceed the checkpoint's completed iteration")
    centroids = summary["refreshed_centroids"]
    resolution = tuple(float(value) for value in np.linalg.norm(image.affine[:3, :3], axis=0))
    template_resolution = tuple(float(value) for value in np.linalg.norm(atlas_topology.image.affine[:3, :3], axis=0))
    options = LesionToadsDescriptorLoopOptions(
        images=images,
        atlas_shape_priors=tuple(atlas_shapes.shape_maps[item.name] for item in priors.structures),
        shape_priors=np.ones(13, dtype=np.float32),
        smoothing_priors=np.ones(13, dtype=np.float32),
        object_types=STANDARD13_TYPES,
        lesion_class=10,
        volume_ratio=float(manifest["settings"]["first_thinning_volume_ratio"]),
        thinning_limit=float(thinning_limit),
        image_resolutions_mm=resolution,
        template_resolutions_mm=template_resolution,
        classification=np.asarray(topology.dataobj),
        template_labels=np.asarray([item.label for item in priors.structures], dtype=np.int16),
        n_iter_max=n_iter_max,
        max_distance=float(manifest["settings"]["max_distance"]),
        algorithm_mode="diffeomorphism",
        relation_limit=3,
        initial_class_centroids=np.asarray(centroids["class_centroids"], dtype=np.float32),
        initial_lesion_centroids=np.asarray(centroids["lesion_centroids"], dtype=np.float32),
        initial_black_hole_centroids=np.asarray(centroids["black_hole_centroids"], dtype=np.float32),
        ventricle_distances=distance_state[..., 1],
        gray_matter_distances=distance_state[..., 3],
        brainstem_distances=distance_state[..., 5],
        atlas_names=tuple(item.name for item in priors.structures),
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
        align_shapes=True,
        init_align_iter=1,
        retain_iteration_details=False,
        initial_distance_maps=LesionInitialDistanceMaps(
            distance_state[..., 0].astype(np.bool_),
            distance_state[..., 1].astype(np.int8),
            distance_state[..., 2].astype(np.bool_),
            distance_state[..., 3].astype(np.int8),
            distance_state[..., 4].astype(np.bool_),
            distance_state[..., 5].astype(np.int8),
        ),
    )
    return run_lesion_toads_descriptor_loop(
        branch, registration, starting_iteration=starting_iteration,
        options=options, max_passes=max_passes,
    )


def persist_plus_lesion_registered_loop_checkpoint(
    registered_branch_manifest: str | Path,
    atlas_descriptor_path: str | Path,
    output_dir: str | Path,
    *,
    n_iter_max: int = 99,
    thinning_limit: float = 3.0,
    source_zero_variance_nan: bool | None = None,
    max_passes: int = 1,
    resume: bool = False,
):
    """Run and persist a bounded continuing loop segment for restart.

    ``max_passes`` creates a checkpoint boundary without declaring source
    convergence. The output is resumable by passing its manifest to this
    function for the following segment.
    """
    from .plus_lesion_registered_branch import PlusLesionRegisteredBranchResult

    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    if isinstance(max_passes, bool) or not isinstance(max_passes, int) or max_passes < 1:
        raise ValueError("max_passes must be a positive integer")
    parent_path = Path(registered_branch_manifest).expanduser().resolve()
    parent, _ = _load_complete_checkpoint(parent_path)
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    if not descriptor.is_file() or descriptor.is_symlink():
        raise FileNotFoundError(f"Missing regular atlas descriptor: {descriptor}")
    settings = {
        "n_iter_max": int(n_iter_max),
        "thinning_limit": float(thinning_limit),
        "source_zero_variance_nan": (
            bool(parent["settings"]["source_zero_variance_nan"])
            if source_zero_variance_nan is None else source_zero_variance_nan
        ),
        "max_passes": int(max_passes),
        "max_distance": float(parent["settings"]["max_distance"]),
        "first_thinning_volume_ratio": float(parent["settings"]["first_thinning_volume_ratio"]),
    }
    signature_payload = {
        "stage": "plus-lesion-registered-loop",
        "parent": {"path": str(parent_path), "sha256": _sha256(parent_path)},
        "atlas": {"path": str(descriptor), "sha256": _sha256(descriptor)},
        "settings": settings,
    }
    signature = hashlib.sha256(
        json.dumps(signature_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    root = Path(output_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "lesion-registered-loop.json"
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Registered-loop checkpoint exists; use resume=True: {manifest_path}")
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old.get("run_signature") != signature:
            raise ValueError("Cannot resume registered loop: inputs or settings differ")
        try:
            cached_manifest, cached_outputs = _load_complete_checkpoint(manifest_path)
            hydrated, _, _ = hydrate_plus_lesion_registered_branch(manifest_path)
            if cached_manifest["stage"] == "plus-lesion-registered-loop":
                return PlusLesionRegisteredBranchResult(
                    manifest_path, cached_outputs["memberships"], cached_outputs["segmentation"],
                    cached_outputs["skeleton"], cached_outputs["availability"],
                    cached_outputs["ordering"], cached_outputs["relation_map"],
                    cached_outputs["summary"], True,
                )
        except (OSError, ValueError, KeyError, RegisteredBranchStateError):
            pass
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no registered-loop manifest exists: {manifest_path}")

    loop = continue_plus_lesion_registered_loop(
        parent_path, descriptor,
        n_iter_max=settings["n_iter_max"],
        thinning_limit=settings["thinning_limit"],
        source_zero_variance_nan=settings["source_zero_variance_nan"],
        max_passes=settings["max_passes"],
    )
    branch = loop.final_branch
    if branch.stop or branch.growing is None or branch.relations is None:
        raise RegisteredBranchStateError(
            "a continuing checkpoint boundary requires growth and refreshed relations; "
            "choose a segment that ends before convergence or the iteration limit"
        )
    parent_summary = json.loads(Path(parent["outputs"]["summary"]["path"]).read_text(encoding="utf-8"))
    recrop_path = Path(parent["inputs"]["registered_recrop_manifest"]["path"]).resolve()
    recrop = json.loads(recrop_path.read_text(encoding="utf-8"))
    modalities_image = load_nifti_ras(Path(recrop["outputs"]["modalities"]["path"]))
    reference = nib.Nifti1Image(
        np.asarray(modalities_image.dataobj)[..., 0],
        modalities_image.affine,
        header=modalities_image.header,
    )
    distance_maps = loop.final_distance_maps
    if distance_maps is None:
        raise RegisteredBranchStateError("GM-enabled loop pass did not refresh distance maps")
    distance_state = np.stack((
        distance_maps.ventricle_mask,
        distance_maps.ventricle,
        distance_maps.gray_matter_mask,
        distance_maps.gray_matter,
        distance_maps.brainstem_mask,
        distance_maps.brainstem,
    ), axis=-1).astype(np.int16, copy=False)
    membership = branch.membership_result.membership
    packed_memberships = np.concatenate((
        membership.class_memberships,
        membership.lesion_membership[..., None],
        membership.black_hole_membership[..., None],
        membership.outlier_membership[..., None],
        membership.lesion_relation[..., None],
    ), axis=-1).astype(np.float32, copy=False)
    diagnostics = np.concatenate((
        membership.denominator[..., None], membership.class_energies,
        membership.lesion_energy[..., None], membership.black_hole_energy[..., None],
        membership.running_objective[..., None],
    ), axis=-1).astype(np.float32, copy=False)
    outputs = {
        "memberships": save_nifti_ras(
            nifti_ras_on_reference(reference, packed_memberships, dtype=np.float32),
            root / "registered_memberships_17ch.nii.gz", dtype=np.float32,
        ),
        "membership_diagnostics": save_nifti_ras(
            nifti_ras_on_reference(reference, diagnostics, dtype=np.float32),
            root / "registered_membership_diagnostics_17ch.nii.gz", dtype=np.float32,
        ),
        "distance_state": save_nifti_ras(
            nifti_ras_on_reference(reference, distance_state, dtype=np.int16),
            root / "initial_distance_state_6ch.nii.gz", dtype=np.int16,
        ),
        "segmentation": save_nifti_ras(
            nifti_ras_on_reference(reference, branch.segmentation, dtype=np.int16),
            root / "registered_segmentation.nii.gz", dtype=np.int16,
        ),
        "skeleton": save_nifti_ras(
            nifti_ras_on_reference(reference, branch.growing.skeleton, dtype=np.int16),
            root / "registered_skeleton.nii.gz", dtype=np.int16,
        ),
        "availability": save_nifti_ras(
            nifti_ras_on_reference(reference, branch.growing.available.astype(np.uint8), dtype=np.uint8),
            root / "registered_availability_13ch.nii.gz", dtype=np.uint8,
        ),
        "ordering": save_nifti_ras(
            nifti_ras_on_reference(reference, branch.growing.ordering, dtype=np.float32),
            root / "registered_ordering.nii.gz", dtype=np.float32,
        ),
        "relation_map": save_nifti_ras(
            nifti_ras_on_reference(reference, branch.relations.relation_map, dtype=np.int16),
            root / "registered_relation_map.nii.gz", dtype=np.int16,
        ),
    }
    ratio = branch.membership_result.energy_ratio
    bounds = branch.registration_handoff.shape_bounds
    completed_iteration = loop.passes[-1].iteration
    summary_output = {
        "spatial_convention": "RAS",
        "class_names": parent_summary["class_names"],
        "transition_objective": float(membership.objective),
        "energy_ratio": float(ratio.energy_ratio),
        "energy_ratio_state": {key: float(getattr(ratio, key)) for key in (
            "previous_energy", "current_energy", "difference", "numerator",
            "denominator", "quotient", "energy_ratio",
        )},
        "refreshed_centroids": {
            "class_centroids": loop.final_class_centroids.tolist(),
            "lesion_centroids": loop.final_lesion_centroids.tolist(),
            "black_hole_centroids": loop.final_black_hole_centroids.tolist(),
        },
        "registration_transform": branch.registration_handoff.transform.tolist(),
        "registration_rotation": branch.registration_handoff.rotation.tolist(),
        "registration_shape_bounds": [{
            "lower_xyz": list(item.lower_xyz), "upper_xyz": list(item.upper_xyz),
            "support_count": item.support_count,
        } for item in bounds],
        "stop": False,
        "stop_reason": None,
        "growing_events": len(branch.growing.events),
        "relation_count": branch.relations.relation_count,
        "relation_diffusion_iterations": branch.relations.diffusion_iterations,
        "relation_membership": branch.relations.relation_membership.astype(np.uint8).tolist(),
        "class_volume": branch.growing.class_volume.tolist(),
        "stage_order": list(loop.passes[-1].stage_order),
        "completed_iteration": int(completed_iteration),
        "source_zero_variance_nan": settings["source_zero_variance_nan"],
        "equivalence_note": "Synthetic checkpoint boundary from selected P39 loop APIs; not complete Java composite equivalence.",
    }
    summary_path = root / "registered_loop_summary.json"
    _json(summary_path, summary_output)
    outputs["summary"] = summary_path
    complete_manifest = {
        "schema_version": 1,
        "stage": "plus-lesion-registered-loop",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "complete",
        "spatial_convention": "RAS",
        "settings": settings,
        "inputs": {
            "parent_checkpoint": {"path": str(parent_path), "sha256": _sha256(parent_path)},
            "registered_recrop_manifest": parent["inputs"]["registered_recrop_manifest"],
            "atlas_descriptor": {"path": str(descriptor), "sha256": _sha256(descriptor)},
        },
        "outputs": {key: {"path": str(path), "sha256": _sha256(path)} for key, path in outputs.items()},
    }
    _json(manifest_path, complete_manifest)
    return PlusLesionRegisteredBranchResult(
        manifest_path, outputs["memberships"], outputs["segmentation"], outputs["skeleton"],
        outputs["availability"], outputs["ordering"], outputs["relation_map"], summary_path, False,
    )
