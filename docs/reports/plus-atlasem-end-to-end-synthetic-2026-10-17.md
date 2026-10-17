# Complete AtlasEM-enabled PLUS synthetic workflow

Date: 2026-10-17  
Test and fixture commit: `ace0cbd6d8fd5e6519d72cde6b16d31fd803f3b3`

## Runnable command

```text
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/pytest -q tests/test_plus_full_workflow.py::test_plus_full_run_with_atlasem_handoff_completes_and_resumes_all_ten_stages
```

The test builds 25³ synthetic PLUS inputs from the committed folded phantom,
uses the committed 24³ two-class AtlasEM atlas in
[`tests/data/synthetic/plus_atlasem_handoff/`](../../tests/data/synthetic/plus_atlasem_handoff/),
then invokes the actual `plus-full-run` CLI. It completes and resumes the
following ten ordered stages:

1. Preprocessing and T1/FLAIR registration
2. AtlasEM and pre-LesionTOADS Dura
3. LesionTOADS crop
4. Initial atlas alignment
5. Registered recrop
6. LesionTOADS iterative workflow
7. LesionTOADS Dura outputs
8. Native-T1 Dura transfer
9. MaSeg cortical separation
10. CRUISE surfaces and thickness outputs

The test verifies all ten stage manifests complete, the parent run is terminal,
AtlasEM/Dura precedes LesionTOADS crop, all emitted NIfTI outputs are finite
RAS `.nii.gz`, and all stages are reused on resume. The initial complete run
and resume took about 10 seconds locally; exact timing depends on the host.
Focused regression covering this test and the existing baseline full-run,
resume, corruption-repair route: **2 passed in 19.21s**.

## Full-suite validation

Test/fixture commit: `ace0cbd6d8fd5e6519d72cde6b16d31fd803f3b3`  
Command:

```text
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-17-plus-atlasem-end2end .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q
```

Result: **1,679 passed, 4 skipped, 0 failed, 30 expected warnings in 774.49s**.
Peak sampled RSS was **14.035 GB**; the 16 GB watchdog did not terminate the
run (`child_exit=0`).

## Evidence boundary

This establishes complete synthetic execution and restart through the optional
AtlasEM/Dura path. The generated atlas is a small integration fixture, not a
replacement for the unavailable selected `ibsr8obj` atlas. It does not validate
human AtlasEM convergence, selected-atlas behavior, anatomy quality, or Java /
MATLAB whole-workflow equivalence. The full PLUS source comparison remains
open.
