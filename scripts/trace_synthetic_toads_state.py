#!/usr/bin/env python3
"""Emit iteration-boundary state summaries for the synthetic registered TOADS loop.

This diagnostic-only trace summarizes the Python implementation. It is intended
to be compared with a corresponding Java trace; it does not establish source
equivalence or validate behavior on real images.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_atlas import (
    load_toads_atlas,
    load_toads_intensity_variance_factors,
    load_toads_optimized_factors,
)
from pymacruise.toads_registered_loop import (
    ToadsRegisteredLoopConfig,
    run_registered_toads_loop,
)
import pymacruise.toads_registered_loop as registered_loop


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/toads_registered_first_iteration_cycle"
ATLAS = ROOT / "tests/data/synthetic/toads_registered_field_refit/atlas-11.txt"


def _float32_bits(value: Any) -> str:
    return f"0x{struct.unpack('<I', np.asarray(value, dtype='<f4').tobytes())[0]:08x}"


def _array_summary(value: Any) -> dict[str, Any]:
    array = np.ascontiguousarray(np.asarray(value))
    summary: dict[str, Any] = {
        "shape": list(array.shape),
        "dtype": str(array.dtype),
        "sha256": hashlib.sha256(array.view(np.uint8)).hexdigest(),
    }
    if array.dtype == np.dtype(np.float32):
        bits = array.view(np.uint32).reshape(-1)
        finite = np.isfinite(array)
        summary["float32_bits_first_12"] = [f"0x{int(bit):08x}" for bit in bits[:12]]
        summary["finite_count"] = int(np.count_nonzero(finite))
        if np.any(finite):
            summary["finite_min"] = float(np.min(array[finite]))
            summary["finite_max"] = float(np.max(array[finite]))
    return summary


def _registered_cycle_summary(
    exact: Any, cycle: Any, index: int, field_snapshot: np.ndarray | None
) -> dict[str, Any]:
    branch = cycle.registered_iteration.branch
    handoff = cycle.registration_control.handoff
    registration = handoff.registration
    growing = branch.growing
    thinning = cycle.next_thinning.thinning if cycle.next_thinning is not None else None
    relations = branch.relations
    final_segmentation = getattr(cycle, "final_segmentation", None)
    if final_segmentation is None and growing is not None:
        final_segmentation = growing.segmentation
    return {
        "iteration": index,
        "stage_order": list(cycle.stage_order),
        "membership": {
            "objective": float(exact.membership_pass.objective),
            "objective_float32_bits": _float32_bits(exact.membership_pass.objective),
            "class_data_energy": _array_summary(exact.membership_pass.class_data_energy),
            "memberships": _array_summary(exact.membership_pass.memberships),
            "denominator": _array_summary(exact.membership_pass.denominator),
            "centroids": _array_summary(exact.centroids),
            "class_standard_deviations": _array_summary(
                getattr(exact, "class_standard_deviations", getattr(exact, "standard_deviations", None))
            ),
            "image_variance_weights": _array_summary(exact.image_variance_weights),
            "image_standard_deviation": _array_summary(exact.image_standard_deviation),
            "centroid_similarity": _array_summary(exact.centroid_similarity),
            "fields": _array_summary(field_snapshot) if field_snapshot is not None else None,
        },
        "registration": {
            "transform": _array_summary(handoff.transform),
            "energy_ratio": float(branch.control.energy_ratio),
            "energy_ratio_float32_bits": _float32_bits(branch.control.energy_ratio),
            "stop": bool(branch.control.stop),
            "run_growing": bool(branch.control.run_growing),
            "update_classification": bool(branch.control.update_classification),
            "refresh_relations": bool(branch.control.refresh_relations),
            "stop_reason": registration.stop_reason if registration is not None else None,
            "iterations": [
                {
                    "transform": _array_summary(item.transform),
                    "energy_float32_bits": _float32_bits(item.energy),
                    "lambda_float32_bits": _float32_bits(item.lambda_value),
                    "relative_change_float32_bits": _float32_bits(item.relative_change),
                    "accepted": bool(item.accepted),
                    "rejected_trials": int(item.rejected_trials),
                }
                for item in (registration.iterations if registration is not None else ())
            ],
        },
        "post_registration": {
            "stage_order": list(branch.stage_order),
            "final_segmentation": (
                _array_summary(final_segmentation) if final_segmentation is not None else None
            ),
            "classification": (
                _array_summary(branch.classification) if branch.classification is not None else None
            ),
            "changed_classification_voxels": branch.changed_classification_voxels,
            "growing": None
            if growing is None
            else {
                "segmentation": _array_summary(growing.segmentation),
                "skeleton": _array_summary(growing.skeleton),
                "ordering": _array_summary(growing.ordering),
                "available": _array_summary(growing.available),
                "class_volumes": _array_summary(growing.class_volumes),
                "event_count": len(growing.events),
                "remaining_queue_count": len(growing.remaining_queue),
            },
            "relations": None
            if relations is None
            else {
                "relation_map": _array_summary(relations.relation.relation_map),
                "relation_membership": _array_summary(relations.relation.relation_membership),
                "boundary_distances": _array_summary(relations.boundary_distances),
                "closest_labels": _array_summary(relations.closest_labels),
                "boundary_diffusion_passes": relations.boundary_diffusion_passes,
                "closest_label_passes": relations.closest_label_passes,
            },
            "next_thinning": None
            if thinning is None
            else {
                "segmentation": _array_summary(thinning.segmentation),
                "skeleton": _array_summary(thinning.skeleton),
                "available": _array_summary(thinning.available),
                "class_volumes": _array_summary(thinning.class_volumes),
                "boundary_counts": _array_summary(thinning.boundary_counts),
                "minimum_class_volumes": _array_summary(thinning.minimum_class_volumes),
                "processed_item_count": len(thinning.processed_items),
            },
        },
    }


def build_trace(
    maximum_iterations: int = 4,
    padding: int = 2,
    arrays_output: Path | None = None,
    relation_limit: int = 3,
) -> dict[str, Any]:
    """Run a deterministic low-cost synthetic trace through the iteration cap."""
    if maximum_iterations < 2:
        raise ValueError("maximum_iterations must be at least 2")
    if padding < 0 or relation_limit < 1:
        raise ValueError("padding must be nonnegative and relation_limit must be positive")
    t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    brain_mask = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    atlas = load_toads_atlas(ATLAS)
    alignment = register_toads_from_atlas_assets(
        t1,
        brain_mask,
        atlas,
        levels=1,
        max_iterations=1,
        padding=padding,
        relation_limit=relation_limit,
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
    )
    category_exponents = load_toads_optimized_factors(ATLAS, [alignment.modality])
    class_indices = {"optimized_csf": 0, "optimized_gm": 1, "optimized_wm": 2}
    exponent_by_class = np.asarray(
        [[
            1.0
            if item.topology.lower() == "mask"
            else category_exponents[0, class_indices[item.topology.lower()]]
            for item in atlas.structures
        ]],
        dtype=np.float32,
    )
    config = ToadsRegisteredLoopConfig(
        category_exponents=category_exponents,
        exponent_by_channel_class=exponent_by_class,
        modality_weights=np.ones(1, dtype=np.float32),
        class_variance_weights=load_toads_intensity_variance_factors(
            atlas.descriptor, [alignment.modality]
        ),
        smoothing_priors=np.ones(len(atlas.structures), dtype=np.float32),
        atlas_shape_priors=np.ones(len(atlas.structures), dtype=np.float32),
        smoothing=0.2 / 6.0 / len(atlas.structures),
        prior_factor=0.5,
        prior_coefficient=2.0 / len(atlas.structures),
        maximum_iterations=maximum_iterations,
        convergence_threshold=0.0,
        relation_limit=relation_limit,
        do_growing=True,
        maximum_polynomial_degree=3,
        maximum_distance=1e-4,
        max_registration_iterations=1,
        registration_energy_tolerance=0.0,
    )
    field_snapshots: dict[int, np.ndarray] = {}
    first_exact_original = registered_loop.continue_registered_unity_startup_through_first_exact_iteration
    following_exact_original = registered_loop.continue_registered_cycle_through_following_exact_iteration

    def capture_first_exact(*args: Any, **kwargs: Any) -> Any:
        exact = first_exact_original(*args, **kwargs)
        field_snapshots[1] = np.asarray(exact.fields, dtype=np.float32).copy()
        return exact

    def capture_following_exact(*args: Any, **kwargs: Any) -> Any:
        exact = following_exact_original(*args, **kwargs)
        field_snapshots[len(field_snapshots) + 1] = np.asarray(kwargs["fields"], dtype=np.float32).copy()
        return exact

    registered_loop.continue_registered_unity_startup_through_first_exact_iteration = capture_first_exact
    registered_loop.continue_registered_cycle_through_following_exact_iteration = capture_following_exact
    try:
        result = run_registered_toads_loop(alignment, config)
    finally:
        registered_loop.continue_registered_unity_startup_through_first_exact_iteration = first_exact_original
        registered_loop.continue_registered_cycle_through_following_exact_iteration = following_exact_original
    startup = result.startup
    if arrays_output is not None:
        arrays_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            arrays_output,
            transformed_shape_priors=np.asarray(alignment.transformed_shape_priors),
            closest_labels=np.asarray(alignment.relation_state.closest_labels),
            segmentation=np.asarray(alignment.segmentation_constructor.segmentation),
            recropped_images_channel_first=np.asarray(
                alignment.template_refresh.recropped_images
            ),
            startup_memberships=np.asarray(
                startup.statistics_membership.membership_pass.memberships
            ),
            first_exact_memberships=np.asarray(
                result.exact_iterations[0].membership_pass.memberships
            ),
        )
    initial_thinning = startup.thinning
    exact_cycles = [
            _registered_cycle_summary(exact, cycle, index, field_snapshots.get(index))
        for index, (exact, cycle) in enumerate(
            zip(result.exact_iterations, result.registered_cycles, strict=True), start=1
        )
    ]
    return {
        "evidence_kind": "python_only_synthetic_diagnostic",
        "source_equivalence_claim": False,
        "real_data_validation": False,
        "fixture": str(FIXTURE.relative_to(ROOT)),
        "shape": list(t1.shape),
        "padding": padding,
        "relation_limit": relation_limit,
        "registration_memberships": {
            "pre_registration": _array_summary(
                alignment.alignment.startup.initial_memberships
            ),
            "post_registration": _array_summary(
                alignment.segmentation_initialization.memberships
            ),
        },
        "membership_inputs": {
            "transformed_shape_priors": _array_summary(
                alignment.transformed_shape_priors
            ),
            "closest_labels": _array_summary(
                alignment.relation_state.closest_labels
            ),
            "segmentation": _array_summary(
                alignment.segmentation_constructor.segmentation
            ),
            "recropped_images_channel_first": _array_summary(
                alignment.template_refresh.recropped_images
            ),
        },
        "iteration_limit": maximum_iterations,
        "convergence_threshold": 0.0,
        "do_growing": config.do_growing,
        "category_exponents": np.asarray(config.category_exponents).tolist(),
        "smoothing": config.smoothing,
        "prior_coefficient": config.prior_coefficient,
        "prior_factor": config.prior_factor,
        "prior_scale": config.prior_scale,
        "stage_order": list(startup.stage_order),
        "initial_state": {
            "fields": _array_summary(startup.fields),
            "objective": float(startup.statistics_membership.membership_pass.objective),
            "objective_float32_bits": _float32_bits(
                startup.statistics_membership.membership_pass.objective
            ),
            "memberships": _array_summary(
                startup.statistics_membership.membership_pass.memberships
            ),
            "centroids": _array_summary(
                alignment.segmentation_initialization.centroids
            ),
            "class_standard_deviations": _array_summary(
                startup.statistics_membership.class_variances.standard_deviations
            ),
            "image_variance_weights": _array_summary(
                startup.statistics_membership.image_variance.image_variance_weights
            ),
            "image_standard_deviation": _array_summary(
                startup.statistics_membership.image_variance.image_standard_deviations
            ),
            "centroid_similarity": _array_summary(
                startup.statistics_membership.centroid_similarity.similarity
            ),
            "segmentation": _array_summary(initial_thinning.segmentation),
            "skeleton": _array_summary(initial_thinning.skeleton),
            "available": _array_summary(initial_thinning.available),
            "class_volumes": _array_summary(initial_thinning.class_volumes),
        },
        "iterations": exact_cycles,
        "terminal": {
            "iteration_count": len(result.exact_iterations),
            "energy_ratio": float(result.terminal_energy_ratio),
            "energy_ratio_float32_bits": _float32_bits(result.terminal_energy_ratio),
            "threshold_reached": result.energy_threshold_reached,
            "iteration_limit_reached": result.iteration_limit_reached,
            "native_classification": (
                _array_summary(result.native_classification)
                if result.native_classification is not None
                else None
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-iterations", type=int, default=4)
    parser.add_argument("--padding", type=int, default=2)
    parser.add_argument("--relation-limit", type=int, default=3)
    parser.add_argument(
        "--arrays-output",
        type=Path,
        default=None,
        help="optionally save membership input arrays as a compressed NPZ diagnostic",
    )
    parser.add_argument("--output", type=Path, default=None, help="write JSON here; stdout if omitted")
    args = parser.parse_args()
    if args.maximum_iterations < 2:
        parser.error("--maximum-iterations must be at least 2")
    if args.padding < 0:
        parser.error("--padding must be non-negative")
    if args.relation_limit < 1:
        parser.error("--relation-limit must be positive")
    payload = json.dumps(
        build_trace(
            args.maximum_iterations,
            padding=args.padding,
            arrays_output=args.arrays_output,
            relation_limit=args.relation_limit,
        ),
        indent=2,
        sort_keys=True,
    ) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
