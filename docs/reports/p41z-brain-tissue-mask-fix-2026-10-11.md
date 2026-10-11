# P41z MaCRUISE+ brain-tissue-mask correction — 2026-10-11

## Source diagnosis

The top-level PLUS workflow and checkpoint composer passed `Atlas_Mask`, which
contains WM labels only, into the MaCRUISE+ tissue-fusion stage. The selected
Java workflow passes `Atlas_BrainTissueMask`, which contains GM, WM, and other
brain-tissue labels. This made every label-100 GM parcel voxel fall outside the
fusion mask, even though the selected MaSeg GM atlas map marked all of them as
GM.

The crop-level selected-source comparison used the same six serialized input
arrays for Java and Python, with the brain-tissue mask cropped from the same
saved crop bounds. Source commit under comparison: `2c28ad12`. The input grid
was RAS, shape 73×91×63.

| Output | Java/Python difference |
|---|---:|
| Fused GM | 0 / 418,509 voxels; bitwise exact |
| Fused WM | 0 / 418,509 voxels; bitwise exact |
| Fused mask | 0 / 418,509 voxels; bitwise exact |

At the 67,677 label-100 voxels, the corrected mask provides support at all
67,677. There are 293 confident-GM and 67,384 fallback voxels, with 0 outside
the fusion mask and 0 confident-WM. Fused GM is nonzero at 67,672 voxels. These
measurements confirm the source-selected fusion behavior on the held-fixed
inputs; they do not establish whole-workflow equivalence.

Machine-readable comparison: [Java/Python fusion result](p41z-brain-mask-java-comparison-2026-10-11.json).

## Corrected iteration-6 checkpoint composition

The full workflow and checkpoint-composition entry points now pass
`Atlas_BrainTissueMask`. A new composed run from saved registered-loop iteration
6 completed its Dura transfer, MaSeg separation, tissue update, and surface
stages in 0.505 GB peak sampled RSS under a 12 GB watchdog. This is a
nonterminal checkpoint composition; the registered-loop energy ratio remains
0.131597 and has not converged.

The corrected output improves GM support from 0% to 0.099%, but still fails the
requested quality checks. GM support between surfaces is 0.099%; WM support
inside the inner surface is 59.45%; the three meshes have four components each
and Euler characteristic −8; median mesh thickness is 0.374 mm. The 5 mm
nominal ribbon bounds and the ≥0.97 tissue-support criterion remain unmet. This
confirms the mask wiring defect was real and that additional independent
quality failures remain. Thresholds were not changed.

The full quality audit is in [the corrected checkpoint report](p41z-mask-corrected-ribbon-quality-2026-10-11.json).

## Validation commands

Focused wiring regression:

```sh
PYTHONPATH=src .venv/bin/pytest -q --tb=short \
  tests/test_plus_full_workflow.py::test_plus_full_workflow_routes_brain_tissue_mask_to_surface_stage
```

The test passed: **1 passed in 0.76 s**. It uses mocked upstream stages and
asserts the full workflow passes the distinct brain-tissue mask path to the
surface workflow.

The focused mask/fusion regression group also passed **6 tests in 0.94 s**:

```sh
PYTHONPATH=src .venv/bin/pytest -q --tb=short \
  tests/test_plus_tissue_update.py tests/test_cortical_separation.py \
  tests/test_plus_full_workflow.py::test_plus_full_workflow_routes_brain_tissue_mask_to_surface_stage
```

Reproduce the held-fixed Java comparison from the repository root:

```sh
scripts/compare_plus_tissue_update_saved_java.sh \
  build/p41z-case-normalization/iteration-0006-surfaces/plus-surfaces \
  build/p41z-case-normalization/brain-mask-java \
  docs/reports/p41z-brain-mask-java-comparison-2026-10-11.json \
  build/p41z-case-normalization/iteration-0006-surfaces/maseg-cortical-separation/Atlas_BrainTissueMask.nii.gz
```

The separate full pytest run in [suite evidence](full-pytest-2026-10-11-lesion-type-normalization.md)
predates this wiring correction. This focused validation does not replace a
new full-suite run after code changes.
