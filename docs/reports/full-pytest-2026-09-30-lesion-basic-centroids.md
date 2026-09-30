# Full pytest evidence: P39zzt generic LesionTOADS centroid modes

- Run: `.venv/bin/python -m pytest -q`
- Result: **618 passed, 2 warnings, in 123.29s**.
- Warnings: two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `8275a1370ab76c8c5beab68518d65ab0372b6771`.
- Focused selected-Java tests: `tests/test_lesion_basic_centroids.py` (**2 passed**).
- Synthetic input and reference provenance: [LesionTOADS variance phantom](../../tests/data/synthetic/lesion_variances/README.md).

The selected generic branch and Python port match raw float32 centroid bits
for `none`, `linked`, `damped`, `prior`, and fallback modes, with and without
field correction. Cases include two channels, a RAYLEIGH class, a class with
no hard-label support, per-channel linked groups, prior centroids, and the
damped-mode empty-class carry-forward. Malformed labels, fields, link maps,
missing prior centroids, and specialized LESION/BLACKHOLE/PowerTable object
types reject explicitly.

This validates the generic centroid update for supplied labels and
memberships. It does not validate specialized centroid branches or the
complete LesionTOADS membership iteration.
