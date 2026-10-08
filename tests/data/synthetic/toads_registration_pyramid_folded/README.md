# Folded-cortex TOADS registration pyramid phantom

This committed, deterministic 48×48×48 RAS fixture contains soft background,
white-matter, folded cortical gray-matter, and CSF fields. The atlas shape
priors are shifted 1.5 voxels from the input memberships so the selected
coarse-to-fine rigid registration takes accepted updates at scales 8, 4, and
2. It complements the smaller four-level startup phantom, which exercises the
degenerate nonfinite-trial edge case.

`tests/test_toads_registration_pyramid_folded.py` calls the actual supplied
VandyCRUISE `MultipleToadDeformableAtlas.updateShapeRegistration()` and
`registerShapesPyramid()` harness and compares every scale's line-search
counts, final six transform float32 values, and last accepted energy/lambda
against Python. This is selected-source synthetic method evidence, not full
TOADS initialization or workflow equivalence.

Regenerate both NIfTIs with:

```bash
./.venv/bin/python tests/data/synthetic/toads_registration_pyramid_folded/generate_fixture.py
```
