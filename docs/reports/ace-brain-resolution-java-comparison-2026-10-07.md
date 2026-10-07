# Brain-resolution MaACE comparison against selected Java — 2026-10-07

## Scope and evidence category

This compares the selected MaCRUISE ACE wrapper with Python on the same
synthetic cropped gray-matter, white-matter, and multi-atlas CSF maps at
`(166, 208, 156)`. The crop is derived from the committed brain-resolution
synthetic base-workflow fixture; no subject MRI is used. It is selected-stage
source-comparison evidence, not end-to-end MaCRUISE equivalence or real-data
validation.

Inputs, selected Java final outputs, Java source hashes, parameter values,
RAS affine, SHA-256 checks, output metrics, and reference-generation command
are recorded in
[`tests/data/synthetic/base_workflow_brain_resolution/ace_maace/provenance.json`](../../tests/data/synthetic/base_workflow_brain_resolution/ace_maace/provenance.json).
The three inputs and three Java `.nii.gz` outputs are committed alongside it.
All images share the RAS geometry and have approximately 0.9375 mm voxel size.
The Java source stays in the existing `tests/java_reference/upstream` corpus.

## Comparison result

The selected Java and Python wrappers produced voxel-identical output arrays
for all 5,386,368 voxels in each output:

| Output | Changed voxels | Maximum absolute error | Mean absolute error | Dice / foreground |
|---|---:|---:|---:|---:|
| Enhanced GM | 0 | 0 | 0 | continuous array equality |
| Cleaned CSF skeleton | 0 | 0 | 0 | Dice 1.0; 74,734 foreground voxels |
| Thinned CSF skeleton | 0 | 0 | 0 | Dice 1.0; 49,378 foreground voxels |

This meets the declared segmentation threshold (Dice ≥0.97) on these outputs.
The arrays were finite. In the separately checked 25-iteration anisotropic
diffusion stage, Java and source-faithful Python each had zero NaNs, there were
zero NaN-mask disagreements, and zero voxels were excluded. The ordinary
Python behavior that replaces zero-tensor divide-by-zero NaNs therefore did
not affect this fixture.

## Source discrepancy found and corrected

The full crop revealed that the selected Java wrapper stores its threshold as
a double, subtracts float WM membership values in double precision, and then
stores the result as float. Python had converted the threshold to float32
before subtraction. That early conversion changed 368,867 of 5,386,368 initial
level-set voxels (maximum absolute difference `1.5258789e-5`, mean
`7.13e-7`) and shifted the thin-skeleton candidate set. Python now performs the
subtraction in float64 before casting the result to float32, matching Java's
operation order. The full-crop output comparison then became exact for all
three outputs.

Focused regression after the correction:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  ./.venv/bin/pytest -q tests/test_ace_reference.py
```

Result: **4 passed in 13.54 s**, sampled peak RSS **0.157 GB**. The separate
full-crop Python-versus-Java execution used the same 25-iteration parameters
and peaked at **1.462 GB** under its 8 GB watchdog. No full repository suite
was run for this localized slice. These focused results cover the working-tree
implementation committed with this report; they are not a full-suite claim.

## Reproduction

Generate the selected Java output from committed inputs with a JDK:

```sh
./.venv/bin/python scripts/generate_java_ace_reference.py \
  --input-root tests/data/synthetic/base_workflow_brain_resolution/ace_maace/inputs \
  --output-dir tests/data/synthetic/base_workflow_brain_resolution/ace_maace/java_reference \
  --place-to-add-csf 0
```

Run the matching Python operation and exact comparison with the same three
inputs and parameters using this command:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_brain_resolution_ace_java.py
```

The command passed: all three final outputs matched voxel-for-voxel, sampled
peak RSS was **1.450 GB**, and the watchdog did not terminate the process. The
script reports hashes and error metrics. Java intermediate stage products are
available in the ignored local build directory from the earlier run, but are
not committed; this corpus keeps only reusable inputs and final goldens.

## Remaining boundary

The evidence validates selected MaACE behavior on a large, analytic synthetic
crop. It does not compare ACE inside the complete base workflow against a
historical MaCRUISE execution, validate every `place_to_add_csf` profile, or
validate human MRI. The corrected double-subtraction behavior is covered by a
focused unit regression; zero-tensor NaN behavior remains explicitly
characterized by the existing 3³ synthetic test.
