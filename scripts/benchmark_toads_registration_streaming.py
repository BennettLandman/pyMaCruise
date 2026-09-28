"""Reproducible small-grid time/temporary-memory probe for TOADS registration."""

from __future__ import annotations

import argparse
import json
import time
import tracemalloc

import numpy as np

from pymacruise.toads_registration_sampling import (
    stream_toads_registration_coefficients,
    stream_toads_registration_energy,
)


def benchmark(size: int) -> dict[str, float | int]:
    dimensions = (size, size, size)
    template = (size + 2, size + 2, size + 2, 3)
    x, y, z, k = np.indices(template, dtype=np.float32)
    atlas = (0.15 + 0.018 * x + 0.014 * y + 0.011 * z + 0.025 * k).astype(np.float32)
    ix, iy, iz = np.indices(dimensions, dtype=np.float32)
    memberships = np.empty(dimensions + (3,), dtype=np.float32)
    memberships[..., 0] = 0.1
    memberships[..., 1] = 0.3 + 0.01 * ix + 0.008 * iy
    memberships[..., 2] = 0.2 + 0.009 * iy + 0.007 * iz
    lower = np.zeros((3, 3), dtype=np.int64)
    upper = np.broadcast_to(np.asarray(dimensions, dtype=np.int64), (3, 3)).copy()
    transform = (0.03, -0.02, 0.025, 0.1, -0.08, 0.04)

    start = time.perf_counter()
    coefficients = stream_toads_registration_coefficients(
        memberships,
        atlas,
        image_resolutions_mm=(1.0, 1.2, 0.8),
        template_resolutions_mm=(1.0, 1.2, 0.8),
        transform=transform,
        class_bounds_low_xyz=lower,
        class_bounds_high_xyz=upper,
    )
    energy = stream_toads_registration_energy(
        memberships,
        atlas,
        image_resolutions_mm=(1.0, 1.2, 0.8),
        template_resolutions_mm=(1.0, 1.2, 0.8),
        transform=transform,
        class_bounds_low_xyz=lower,
        class_bounds_high_xyz=upper,
    )
    elapsed = time.perf_counter() - start

    # Measure Python-visible temporary allocations in a separate call so
    # tracing overhead does not distort the elapsed-time measurement above.
    tracemalloc.start()
    stream_toads_registration_coefficients(
        memberships,
        atlas,
        image_resolutions_mm=(1.0, 1.2, 0.8),
        template_resolutions_mm=(1.0, 1.2, 0.8),
        transform=transform,
        class_bounds_low_xyz=lower,
        class_bounds_high_xyz=upper,
    )
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    if not np.isfinite(energy) or not np.isfinite(coefficients.cost):
        raise RuntimeError(f"synthetic registration benchmark failed for size {size}")
    return {
        "size_xyz": size,
        "classes": 3,
        "voxels": size**3,
        "coefficient_plus_energy_seconds": elapsed,
        "python_tracemalloc_peak_mib": peak / (1024.0**2),
        "energy": float(energy),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=(6, 8, 10))
    args = parser.parse_args()
    if any(size < 2 for size in args.sizes):
        parser.error("all sizes must be at least 2")
    report = {
        "purpose": "synthetic-only streaming TOADS sampler cost/memory probe",
        "measurements": [benchmark(size) for size in args.sizes],
        "memory_note": "tracemalloc measures Python-visible allocations during the sampler call; input arrays are allocated before measurement",
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
