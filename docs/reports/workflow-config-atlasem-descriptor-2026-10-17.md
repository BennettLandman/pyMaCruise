# Config-driven PLUS AtlasEM descriptor handoff

Date: 2026-10-17  
Implementation commit: `3a89195ce6ee71c387d7841bc0d0d765c009d9cc`

## Change

The version-1 workflow JSON can now set optional `atlasem_descriptor_path` for
PLUS. Config-relative paths resolve beside the JSON file. Base configs reject
the field. `workflow-preflight` parses the selected AtlasEM descriptor, checks
that every topology/shape resource exists as `.nii.gz`, verifies 3D RAS headers,
checks descriptor dimensions and common affine grid, and records descriptor
and payload SHA-256 values without materializing the image voxel arrays. The
initial manifest signature includes these resource hashes. `workflow-run`
passes the selected descriptor to the restartable PLUS runner, activating the
AtlasEM→Dura stage before LesionTOADS cropping.

The bundled LesionTOADS descriptor remains the pipeline atlas and the bundled
base TOADS descriptor remains unchanged. A generated descriptor in tests
validates config plumbing and resource contracts; it is not a replacement for
the selected `ibsr8obj` resource and does not validate selected-atlas behavior.

## Focused validation

Command:

```text
PYTHONPATH=.:src PYTHONNOUSERSITE=1 .venv/bin/pytest -q tests/test_workflow_preflight_cli.py
```

Result: **16 passed, 0 failed, 0 skipped in 3.10s**. The cases include relative
descriptor resolution, PLUS runner forwarding, base-profile rejection, missing
payload rejection, mismatched atlas-grid rejection, hashed manifest resources,
and existing config/preflight behavior. Expected orientation warnings came
from the synthetic LAS T1/MaSeg/mask fixtures.

## Full-suite validation

Implementation commit tested: `3a89195ce6ee71c387d7841bc0d0d765c009d9cc`  
Command:

```text
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-17-workflow-atlasem-config .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q
```

Result: **1,678 passed, 4 skipped, 0 failed, 30 expected warnings in 699.17s**.
Peak sampled RSS was **13.843 GB**; watchdog limit was 16 GB and it did not
terminate the run (`child_exit=0`).

## Evidence boundary

This validates the configuration interface, preflight resource checks, and
runner dispatch using synthetic AtlasEM resources. It does not validate the
selected `ibsr8obj` atlas or source-selected AtlasEM/PLUS numerical equivalence.
