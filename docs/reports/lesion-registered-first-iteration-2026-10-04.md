# Registered LesionTOADS first loop update

- Implementation: `run_lesion_toads_registered_first_iteration`
- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_first_iteration.py tests/test_lesion13_finite_grow_thin.py`
- Focused result: **2 passed, 0 failed, 0 skipped in 91.07s**.
- Tested implementation commit: `150e68df478d6248874b596d197f0c2ddfab4c7c`.

The composition starts with P39bm's first thinning result and runs the
selected `computeCentroids` → `computeVariances` → `computeChannelWeights` →
GM-enabled standard-13 `computeMemberships` → float32 energy-ratio sequence.
It carries the prior maps, shape state, relations, and distances through those
operations. Source zero-variance mode is used whenever refreshed deviations
contain zero. The synthetic test directly recomposes the existing components
and checks centroid/statistic/weight/membership/objective/ratio equality and
stage order. The companion P39bm and existing 13-class test cover selected
Java source-excerpt/harness methods individually.

This fixture supplies unit tissue exponents and unit shape/smoothing weights;
those are synthetic seam inputs, not a new claim that the bundled descriptor's
optimized factors have been mapped. The implementation adds **0 selected
TOADS Java source lines** to the coarse source-method accounting. It does not
include field fitting, atlas registration, convergence branching, growth,
whole LesionTOADS equivalence, or real-image validation.

Full-suite evidence is in the [P39bn regression report](full-pytest-2026-10-04-lesion-registered-first-iteration.md).
