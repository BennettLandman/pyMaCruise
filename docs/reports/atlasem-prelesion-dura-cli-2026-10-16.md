# AtlasEM pre-LesionTOADS Dura command-line adapter — 2026-10-16

Added `pymacruise atlasem-prelesion-dura` for the completed AtlasEM-output → Dura handoff. It accepts an existing class-last membership NIfTI, matching T1 and binary brain mask, JSON class names in channel order, iteration limit, and output directory. It writes the three-channel Dura input, AtlasEM classification, and seven Dura products as `.nii.gz`, plus optional RAS provenance JSON.

The command help states that it consumes already-computed memberships. It does not parse the AtlasEM descriptor, load the selected `ibsr8obj` assets, estimate registration, or run the AtlasEM EM loop; it must not be presented as a complete AtlasEM command.

Focused command on implementation commit `b634c7cf441d10920fba55369b600e8be339bbc8`:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_cli.py tests/test_plus_atlas_em_dura_integration.py
```

Result: **2 passed, 0 failed, 0 skipped in 5.26 s**, with three expected warnings for reorienting LAS T1, memberships, and brain mask to RAS. The synthetic CLI test verifies all nine outputs, compressed NIfTI paths, output finiteness, RAS affines, integer segmentation dtype, and provenance settings. This is CLI and synthetic handoff evidence; it is not AtlasEM optimizer parity or full PLUS source equivalence. No full-suite run was performed.
