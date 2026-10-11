# Human MaACE held-fixed Java comparison — 2026-10-11

## Scope and result

This compares the selected Java MaCRUISE ACE wrapper with the already-saved
ACE outputs from the corrected, boundary-enabled human base tail. It reuses
only the saved fused GM, fused WM-fill, cleaned MaSeg CSF boundary, and Python
ACE products; it does not rerun imaging preparation, TOADS, the base workflow,
or any surface stage. Inputs and outputs are on the same 141×173×131 RAS crop
at 1 mm spacing, with affine translation (−70.5, −85.5, −16.5) mm.

The wrapper settings match the saved human run: 25 diffusion iterations,
threshold 126.9, Gaussian alpha 1.8, tensor alpha 1.6, and mode 2 (merge the
supplied boundary after diffusion). Java was run as a direct child of the
16 GB RSS watchdog with `-Xmx4g`; peak sampled Java RSS was **0.663 GB** and it
exited successfully. The selected Java sources and all three fixed inputs are
hashed in the [machine-readable report](../../build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/comparison.json).

| ACE product | Different voxels | MAE | Maximum absolute error | Binary Dice |
|---|---:|---:|---:|---:|
| Enhanced GM | 2 / 3,195,483 | 0.00015154 | 247.4375 | — |
| Skeleton | 3 / 3,195,483 | 0.00023846 | 254 | 0.9999899 |
| Thinned skeleton | 11,165 / 3,195,483 | 0.8874746 | 254 | 0.9411882 |

The two enhanced-GM mismatches are at crop voxels `(51,34,71)` (Java 0,
Python 247.4375) and `(110,115,90)` (Java 236.79486, Python 0). The skeleton
mismatches occur at those locations and `(52,158,80)`. Thus the enhanced-GM
aggregate MAE is small but conceals two large isolated differences; the
thinned skeleton shows a broader topology difference and falls below the
user's 0.97 Dice threshold if assessed as a binary mask. This is a concrete
held-fixed source-comparison failure signal, not a diagnosis of the first
algorithmic stage that diverges.

## Interpretation and limits

The human surface-extraction comparison already shows that Python and selected
Java produce the same 143 connected components from the same saved level sets.
This ACE replay instead finds a small upstream skeleton mismatch that expands
to 11,165 voxels after skeleton thinning. It is a plausible contributor to the
poor human ribbon, but these measurements do not prove causation. The selected
TGDM comparison on the same human run passed its held-fixed surface gates; this
ACE comparison does not establish full-workflow equivalence, identify the
first divergent internal ACE field, validate clinical performance, or provide
a historical whole-pipeline reference.

The next useful bounded diagnostic is a stage-by-stage ACE trace comparison
using Java's saved intermediates and the corresponding Python `ace.py` stages,
with special attention to the two enhanced-GM/skeleton coordinates and to
`ThinACE`'s sensitivity to the three seed discrepancies. Keep the current
saved upstream maps fixed while localizing that difference.

## Reproduction

From the PyMaCRUISE repository root, the exact saved inputs, compiled harness,
and binary output are under
`build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/`.
After compiling `tests/java_reference/ACEReferenceHarness.java` against
`tests/java_reference/stubs` and `tests/java_reference/upstream`, the direct
input preparation is:

```sh
PYTHONPATH=src .venv/bin/python \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/prepare.py
```

Then run Java directly under the watchdog:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms256m -Xmx4g \
  -cp build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/classes \
  ACEReferenceHarness \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/ace_input.bin \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/ace_output.bin \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/intermediate \
  2
```

After Java exits, the held-fixed arrays and summary report are regenerated
without rerunning Java by:

```sh
PYTHONPATH=src .venv/bin/python \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/compare_saved.py
```

This is real-human-data held-fixed evidence for ACE only. It is not a full
workflow rerun, a claim of whole MaCRUISE equivalence, or clinical validation.
