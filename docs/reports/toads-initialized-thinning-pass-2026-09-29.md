# P18az — initialized TOADS thinning pass

## Result

Added toads_propagate_thinning to compose the selected exact-membership
source path: initialize skeleton, availability, per-class volumes and
boundary priorities; compute minVolume[k] = max(100, volumeRatio *
volume[k]) with float32 arithmetic; and run the queue to the source
propagateThinning loop stop condition.

The synthetic reference is a deterministic 12×12×12 grid with four label
classes and 1,000 interior voxels. All four classes begin above 100 voxels.
With volumeRatio=0.8, class 2 uses the 100 floor, while the other three
minimums use the ratio term. Every one of the 1,000 initialized boundary items
is popped; the final heap is empty. The Python pass matches the Java reference
exactly for:

- initial class-volume float32 bits, minimum-volume bits and boundary counts
  [250, 250, 250, 250];
- all 1,000 popped priority bits and (x, y, z, class) payloads;
- final segmentation, skeleton and availability arrays;
- final class-volume float32 bits and residual queue.

The endpoint-counting quirk already identified in P18av remains preserved;
the selected source counter includes all 26 non-center neighbors, making the
Nngb == 1 branch unreachable on these interior neighborhoods.

## Validation evidence

| Check | Result |
|---|---|
| Java source/LUT/heap reference | 1,000 queue items; exact discrete arrays, queue sequence and class-volume bits |
| Focused test | 1 passed in 0.62 s |
| Full pytest | 361 passed in 74.01 s |
| Tested implementation commit | 2d5418aa32584ad2de5e6c9dc9df9085ef6265ba |
| Full-suite command | .venv/bin/pytest -q from pymacruise/ |
| Environment | pytest 9.1.1, Python 3.13.2 |
| Run date | 2026-09-29 |
| Java harness SHA-256 | bd49c7e8849cae45f8088f4d0b898cb20e4aa5a2f0653f5a9cfaa237b8452123 |
| Fixture SHA-256 | 439a21b35a4f71912afe7ae2a766d7c59972f5cfc537acf6175739ed2d1d99d9 |

The committed generator is
[scripts/generate_java_toads_initialized_thinning_pass_reference.sh](../../scripts/generate_java_toads_initialized_thinning_pass_reference.sh);
the comparison test and fixture are
[tests/test_toads_initialized_thinning_pass.py](../../tests/test_toads_initialized_thinning_pass.py)
and
[tests/data/synthetic/toads_initialized_thinning_pass/java_reference.json](../../tests/data/synthetic/toads_initialized_thinning_pass/java_reference.json).

## Scope limits

This is a complete initialized thinning pass for one deterministic synthetic
grid and selected exact-membership settings. It is not evidence for all
possible topology configurations, growing, deformable atlas registration,
Dura integration, or complete TOADS segmentation. No actual imaging data was
used. P18 and P34 remain open.
