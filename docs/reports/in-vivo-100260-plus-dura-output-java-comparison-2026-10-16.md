# Human PLUS LesionTOADS-to-Dura output comparison — 2026-10-16

## Evidence category and scope

This is a selected-source, held-fixed **LesionTOADS output-stage comparison**
on the human-derived case-100260 registered crop. The selected Java
`LesionToadSegmentation.generateDuraRemovalOutputs` method received the exact
Python-generated memberships and segmentation; its four outputs were then
compared with the saved workflow products after the wrapper's ×255 scale.
This does not compare the complete LesionTOADS iteration or an independently
composed original PLUS workflow.

## Inputs and source

- Workflow: ignored `build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/`.
- Exact grid: 164×207×196 on RAS; all maps share the registered membership
  affine. The 17-channel membership file supplies 13 class maps, lesion
  membership (channel 13), and black-hole membership (channel 14); the
  accompanying integer segmentation supplies selected class indices.
- Membership SHA-256:
  `b60c1d7af18921fcfedab526aa45a6d89c00c1b7fec6ce3b216bbb76d1a094a2`.
- Segmentation SHA-256:
  `3e79bf307eec45ec2b4fe410d2ba527160b60b787b04979e56912509d18c300b`.
- Selected source file:
  `VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
  SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
  The test harness is a minimal copy of the selected method body, not a launch
  of the MIPAV/JIST application.

## Results

Each channel has **6,653,808** voxels. After Java's selected method and the
wrapper's float32 ×255 scaling, all four outputs match their saved Python
workflow NIfTIs bit-for-bit:

| Output | Differing voxels | Mean / max absolute error | Nonzero voxels |
|---|---:|---:|---:|
| CSF | 0 | 0 / 0 | 2,413,069 |
| GM | 0 | 0 / 0 | 2,248,345 |
| WM fill | 0 | 0 / 0 | 1,935,772 |
| WM mask | 0 | 0 / 0 | 749,585 |

The selected-source replay peaked at **1.241 GB RSS** under an 8 GB watchdog.
No human images were copied into the repository. Exact machine metrics and
all input/output hashes are in the [comparison JSON](in-vivo-100260-plus-dura-output-java-comparison-2026-10-16.json).

## Reproduction

From the repository root:

```bash
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_human_plus_dura_outputs_java.py \
  --workflow-root build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14 \
  --output-dir build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/java-dura-output-comparison \
  --report build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/java-dura-output-comparison/comparison.json
```

The existing focused synthetic/reference test `tests/test_lesion_dura_outputs.py`
also passes: **3 passed in 3.20 s**. The latest full suite on the unchanged
production implementation passed **1,582 passed, 4 skipped, 0 failed** on
tested tree commit `857370c002164a5775633117059fe35ecdd7aaa8`; see the [suite
report](full-pytest-2026-10-16-ace-human-tgdm.md).
