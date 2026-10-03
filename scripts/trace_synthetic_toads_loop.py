#!/usr/bin/env python3
"""Trace the selected registered TOADS loop on its committed 13-cubed phantom.

This is a diagnostic for integration/runtime characterization, not a new
algorithm implementation or source-equivalence claim. It times existing
startup, exact-membership, and registered-branch functions in place.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import replace
import hashlib
from pathlib import Path
import time
from typing import Any, Callable, Iterator

import nibabel as nib
import numpy as np

from pymacruise.image_io import nifti_ras_on_reference, normalize_image_to_ras
from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_atlas import (
    load_toads_atlas,
    load_toads_intensity_variance_factors,
    load_toads_optimized_factors,
)
from pymacruise.toads_shape_volume import transform_toads_shape_priors
import pymacruise.toads_registered_loop as registered_loop


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/toads_registered_first_iteration_cycle"
ATLAS = ROOT / "tests/data/synthetic/toads_registered_field_refit/atlas-11.txt"


def _timed(name: str, timings: dict[str, list[float]]) -> Callable[..., Any]:
    def decorate(function: Callable[..., Any]) -> Callable[..., Any]:
        def wrapped(*args: Any, **kwargs: Any) -> Any:
            start = time.perf_counter()
            try:
                return function(*args, **kwargs)
            finally:
                timings.setdefault(name, []).append(time.perf_counter() - start)

        return wrapped

    return decorate


@contextmanager
def _measure_loop_stages(timings: dict[str, list[float]]) -> Iterator[None]:
    names = (
        "compute_registered_toads_unity_field_startup",
        "continue_registered_unity_startup_through_first_exact_iteration",
        "continue_registered_first_exact_iteration_through_next_thinning",
        "continue_registered_cycle_through_following_exact_iteration",
        "continue_registered_membership_through_branch",
    )
    originals = {name: getattr(registered_loop, name) for name in names}
    try:
        for name in names:
            setattr(registered_loop, name, _timed(name, timings)(originals[name]))
        yield
    finally:
        for name, function in originals.items():
            setattr(registered_loop, name, function)


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "build/synthetic-toads-trace.nii.gz")
    parser.add_argument(
        "--fixture", type=Path, default=FIXTURE,
        help="directory containing synthetic t1.nii.gz and brain-mask.nii.gz",
    )
    parser.add_argument(
        "--atlas", type=Path, default=ATLAS,
        help="source-format atlas descriptor to use for the registered loop",
    )
    parser.add_argument("--maximum-iterations", type=_positive_int, default=10)
    parser.add_argument("--convergence-threshold", type=float, default=1e-4)
    parser.add_argument("--registration-levels", type=_positive_int, default=1)
    parser.add_argument("--registration-iterations", type=int, default=1)
    parser.add_argument("--loop-registration-iterations", type=_positive_int, default=1)
    parser.add_argument("--padding", type=int, default=2)
    parser.add_argument("--relation-limit", type=_positive_int, default=3)
    parser.add_argument("--maximum-distance", type=float, default=1e-4)
    parser.add_argument(
        "--intensity-variant",
        choices=("original", "constant-volume", "constant-mask", "linear-x", "linear-xyz", "original-scaled"),
        default="original",
    )
    parser.add_argument(
        "--shape-variant",
        choices=("descriptor", "uniform-one", "topology-one-hot"),
        default="descriptor",
        help="loop-time shape-prior probe after the ordinary descriptor-based alignment",
    )
    args = parser.parse_args()

    if args.registration_iterations < 0:
        parser.error("--registration-iterations must be nonnegative")
    if not np.isfinite(args.convergence_threshold) or args.convergence_threshold < 0:
        parser.error("--convergence-threshold must be finite and nonnegative")
    if not np.isfinite(args.maximum_distance) or args.maximum_distance < 0:
        parser.error("--maximum-distance must be finite and nonnegative")
    output = args.output.expanduser().resolve()
    if not output.name.lower().endswith(".nii.gz"):
        parser.error("--output must end with .nii.gz")
    if output.exists():
        parser.error(f"refusing to overwrite existing output: {output}")

    fixture = args.fixture.expanduser().resolve()
    atlas_path = args.atlas.expanduser().resolve()
    t1 = nib.load(str(fixture / "t1.nii.gz"))
    brain_mask = nib.load(str(fixture / "brain-mask.nii.gz"))
    atlas = load_toads_atlas(atlas_path)
    original_values = np.asarray(t1.dataobj, dtype=np.float32)
    mask_values = np.asanyarray(brain_mask.dataobj) > 0
    coordinates = np.indices(t1.shape, dtype=np.float32)
    intensity_variants = {
        "original": original_values,
        "constant-volume": np.full(t1.shape, 1000.0, dtype=np.float32),
        "constant-mask": np.where(mask_values, 1000.0, 0.0).astype(np.float32),
        "linear-x": np.float32(1000.0) + np.float32(20.0) * coordinates[0],
        "linear-xyz": (
            np.float32(1000.0)
            + np.float32(20.0) * coordinates[0]
            + np.float32(10.0) * coordinates[1]
            + np.float32(5.0) * coordinates[2]
        ),
        "original-scaled": (original_values * np.float32(0.001)).astype(np.float32),
    }
    if args.intensity_variant != "original":
        t1 = nib.Nifti1Image(intensity_variants[args.intensity_variant], t1.affine, t1.header.copy())
    start = time.perf_counter()
    alignment = register_toads_from_atlas_assets(
        t1,
        brain_mask,
        atlas,
        levels=args.registration_levels,
        max_iterations=args.registration_iterations,
        padding=args.padding,
        relation_limit=args.relation_limit,
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
    )
    if args.shape_variant != "descriptor":
        topology = alignment.atlas_data.topology
        if args.shape_variant == "uniform-one":
            maps = {
                structure.name: np.ones_like(topology, dtype=np.float32)
                for structure in atlas.structures
            }
        else:
            maps = {
                structure.name: (topology == structure.label).astype(np.float32)
                for structure in atlas.structures
            }
        atlas_data = replace(alignment.atlas_data, shapes=maps)
        shapes = np.stack(tuple(maps[structure.name] for structure in atlas.structures), axis=-1)
        transformed = transform_toads_shape_priors(
            shapes,
            image_dimensions_xyz=tuple(int(value) for value in alignment.transformed_shape_priors.shape[:3]),
            image_resolutions_mm=alignment.alignment.setup.signed_image_resolutions_mm,
            template_resolutions_mm=alignment.alignment.setup.template_resolutions_mm,
            transform=alignment.post_refresh_registration.final_transform,
            class_bounds=alignment.post_refresh_registration.final_shape_bounds,
        )
        alignment = replace(
            alignment,
            atlas_data=atlas_data,
            transformed_shape_priors=transformed,
        )
    alignment_seconds = time.perf_counter() - start

    category_exponents = load_toads_optimized_factors(atlas_path, [alignment.modality])
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
    config = registered_loop.ToadsRegisteredLoopConfig(
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
        maximum_iterations=args.maximum_iterations,
        convergence_threshold=args.convergence_threshold,
        relation_limit=args.relation_limit,
        do_growing=True,
        maximum_polynomial_degree=3,
        maximum_distance=args.maximum_distance,
        max_registration_iterations=args.loop_registration_iterations,
    )
    stage_timings: dict[str, list[float]] = {}
    loop_start = time.perf_counter()
    with _measure_loop_stages(stage_timings):
        result = registered_loop.run_registered_toads_loop(alignment, config)
    loop_seconds = time.perf_counter() - loop_start

    if result.native_classification is None:
        raise RuntimeError("selected registered loop returned without terminal classification")
    ras_t1 = normalize_image_to_ras(t1, source_name="synthetic N4-corrected T1 reference")
    output_image = nifti_ras_on_reference(ras_t1, result.native_classification, dtype=np.uint8)
    output.parent.mkdir(parents=True, exist_ok=True)
    nib.save(output_image, str(output))
    digest = hashlib.sha256(output.read_bytes()).hexdigest()

    print(f"fixture={fixture} shape={t1.shape} atlas={atlas_path}")
    print(
        "settings="
        f"maximum_iterations:{args.maximum_iterations},threshold:{args.convergence_threshold:g},"
        f"registration_levels:{args.registration_levels},registration_iterations:{args.registration_iterations},"
        f"loop_registration_iterations:{args.loop_registration_iterations},relation_limit:{args.relation_limit}"
    )
    print(f"probe.intensity_variant={args.intensity_variant} shape_variant={args.shape_variant}")
    print(
        "source_parameters="
        f"category_exponents:{category_exponents.tolist()},prior_factor:{config.prior_factor:g},"
        f"prior_scale:{config.prior_scale:g},prior_coefficient:{config.prior_coefficient:g},"
        f"smoothing:{config.smoothing:.9g},relation_limit:{config.relation_limit}"
    )
    print(f"elapsed.alignment_seconds={alignment_seconds:.3f}")
    print(f"elapsed.loop_startup_seconds={sum(stage_timings.get('compute_registered_toads_unity_field_startup', [])):.3f}")
    print(f"elapsed.loop_total_seconds={loop_seconds:.3f}")
    for index, (exact, cycle) in enumerate(
        zip(result.exact_iterations, result.registered_cycles, strict=True), start=1
    ):
        branch = cycle.registered_iteration.branch
        if index == 1:
            membership_seconds = sum(
                stage_timings.get("continue_registered_unity_startup_through_first_exact_iteration", [])
            )
            branch_seconds = sum(
                stage_timings.get("continue_registered_first_exact_iteration_through_next_thinning", [])
            )
        else:
            membership_seconds = stage_timings[
                "continue_registered_cycle_through_following_exact_iteration"
            ][index - 2]
            branch_seconds = stage_timings[
                "continue_registered_membership_through_branch"
            ][index - 2]
        registration = cycle.registration_control.handoff.registration
        registration_iterations = len(registration.iterations) if registration is not None else 0
        print(
            f"iteration.{index}=objective:{float(exact.membership_pass.objective):.9g},"
            f"energy_ratio:{float(branch.control.energy_ratio):.9g},"
            f"membership_update_seconds:{membership_seconds:.3f},"
            f"registration_branch_seconds:{branch_seconds:.3f},"
            f"registration_steps:{registration_iterations},stop:{str(branch.control.stop).lower()}"
        )
    print(
        f"terminal=iterations:{len(result.exact_iterations)},"
        f"energy_ratio:{float(result.terminal_energy_ratio):.9g},"
        f"threshold_reached:{str(result.energy_threshold_reached).lower()},"
        f"iteration_limit_reached:{str(result.iteration_limit_reached).lower()}"
    )
    data = np.asanyarray(nib.load(str(output)).dataobj)
    print(
        f"output={output} sha256={digest} shape={data.shape} dtype={data.dtype} "
        f"orientation={nib.aff2axcodes(nib.load(str(output)).affine)} labels={np.unique(data).tolist()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
