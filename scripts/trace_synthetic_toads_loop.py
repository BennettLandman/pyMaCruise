#!/usr/bin/env python3
"""Trace the selected registered TOADS loop on its committed 13-cubed phantom.

This is a diagnostic for integration/runtime characterization, not a new
algorithm implementation or source-equivalence claim. It times existing
startup, exact-membership, and registered-branch functions in place.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
from pathlib import Path
import time
from typing import Any, Callable, Iterator

import nibabel as nib
import numpy as np

from pymacruise.image_io import nifti_ras_on_reference, normalize_image_to_ras
from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_atlas import load_toads_atlas
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
    parser.add_argument("--maximum-iterations", type=_positive_int, default=10)
    parser.add_argument("--convergence-threshold", type=float, default=1e-4)
    parser.add_argument("--registration-levels", type=_positive_int, default=1)
    parser.add_argument("--registration-iterations", type=int, default=1)
    parser.add_argument("--loop-registration-iterations", type=_positive_int, default=1)
    parser.add_argument("--padding", type=int, default=2)
    parser.add_argument("--relation-limit", type=_positive_int, default=1)
    parser.add_argument("--maximum-distance", type=float, default=1e-4)
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

    t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    brain_mask = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    atlas = load_toads_atlas(ATLAS)
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
    alignment_seconds = time.perf_counter() - start

    category_exponents = np.asarray([[2.0, 2.0, 2.0]], dtype=np.float32)
    exponent_by_class = np.asarray(
        [[1.0 if item.topology.lower() == "mask" else 2.0 for item in atlas.structures]],
        dtype=np.float32,
    )
    config = registered_loop.ToadsRegisteredLoopConfig(
        category_exponents=category_exponents,
        exponent_by_channel_class=exponent_by_class,
        modality_weights=np.ones(1, dtype=np.float32),
        class_variance_weights=alignment.initial_variances.standard_deviations,
        smoothing_priors=np.ones(len(atlas.structures), dtype=np.float32),
        atlas_shape_priors=np.ones(len(atlas.structures), dtype=np.float32),
        smoothing=0.2 / 6.0 / len(atlas.structures),
        prior_factor=0.2,
        prior_coefficient=2.0 / len(atlas.structures),
        maximum_iterations=args.maximum_iterations,
        convergence_threshold=args.convergence_threshold,
        relation_limit=args.relation_limit,
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

    print(f"fixture={FIXTURE.relative_to(ROOT)} shape={t1.shape} atlas={ATLAS.relative_to(ROOT)}")
    print(
        "settings="
        f"maximum_iterations:{args.maximum_iterations},threshold:{args.convergence_threshold:g},"
        f"registration_levels:{args.registration_levels},registration_iterations:{args.registration_iterations},"
        f"loop_registration_iterations:{args.loop_registration_iterations},relation_limit:{args.relation_limit}"
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
