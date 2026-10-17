# Initial external BASE test handoff — 2026-10-16

Historical record: the initial v3.1.0 reference below was superseded by the
owner-selected v3.2.0 container. See the [v3.2.0 inspection record](external-v320-reference-inspection-2026-10-16.md)
and current walkthrough for the revised reference and command.

The [walkthrough](../testing/full-base-walkthrough.md) includes Linux/Docker
prerequisites, venv creation, package installation, input hashes, both full
workflow launch commands, Python resume, final scientific-output comparison,
and return-file instructions for a read-only GitHub user.

## External image metadata verification

The public Docker registry returned the following for
`masidocker/spiders:MaCRUISE_v3_1_0`:

- Manifest digest:
  `sha256:9f04f42ea7b2756ec2e38da8ab04992053079393f46a34cc739b596e70097d9f`
- Manifest type: Docker distribution v2.
- Configuration platform: `linux/amd64`.
- Entrypoint: null; default command: `/bin/bash`.
- Sum of compressed layer sizes: approximately 9.483 GB.

Only the manifest and configuration metadata were fetched. The image was not
downloaded or executed. The launch command is based on the public MASILab
MaCRUISE README, using the immutable digest above. It is the complete external
BASE entry point, rather than the adapted selected-mask MATLAB driver.

## Portable comparator

`scripts/compare_external_base_run.py` discovers the native external output
tree and uses the Python manifest's active stage paths. It verifies Python
output hashes and original input identity, enforces an explicit native
volume/mesh coordinate contract, and reports final segmentation Dice, native
level-set errors, three symmetric triangle-surface distance metrics, and
volume/surface thickness summaries. Missing or incompatible products produce
an error report. Numeric label/surface failures remain visible.

The nine synthetic tests cover active-versus-stale attempts, identical outputs,
a known surface displacement, oblique axis-only mapping, physical-grid
rejection, incorrect raw input, output corruption, missing labels, failed
reference execution, and a missing-product error report. These cases are
synthetic comparator tests, not execution of the original imaging application.

Focused command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q tests/test_external_base_comparison.py \
  tests/test_validation.py tests/test_source_run_manifest.py \
  tests/test_source_reference_base_driver.py tests/test_source_reference_plus_driver.py \
  tests/test_source_reference_mask_helper.py
```

Result: **24 passed**, including the nine comparator tests. Documentation's
local links and Bash/Python syntax were also checked.

## Unverified scope

Fresh Linux dependency installation, complete container execution, the exact
export layout/frame contract of this immutable runtime, detailed independent
intermediate mapping, and vertexwise thickness correspondence remain
unverified. External skull stripping and the Python supplied-mask interface
can differ. Scalar/thickness acceptance limits are not universally defined.
The comparator therefore never reports whole-port acceptance, even if its
segmentation and surface gates pass. PLUS's complete selected source stage
chain remains a separate open implementation/asset issue.
