# Synthetic PLUS AtlasEM handoff fixture

This generated two-class atlas exists only to exercise the complete optional
AtlasEM→Dura handoff in the restartable PLUS integration and resume path. It is
not the selected `ibsr8obj` AtlasEM resource, does not represent realistic
brain anatomy, and cannot establish Java/MATLAB source equivalence.

The 24×24×24 topology and shape priors are compressed RAS NIfTI files on a
1-mm grid. The size supports the configured four-level registration pyramid.
The expected CLI test creates synthetic 25³ T1, FLAIR, brain mask, and MaSeg
inputs from the committed `plus_pre_tgdm_composition` folded phantom and runs
`plus-full-run` through all ten stages, then resumes and verifies all stages.

Regenerate descriptor and payloads with:

```text
.venv/bin/python tests/data/synthetic/plus_atlasem_handoff/generate_fixture.py
```

`SHA256SUMS` pins the generated atlas files. `provenance.json` identifies the
fixture as synthetic-only.
