# Integrated second-pass TOADS comparison against Java

## Result

On the corrected seven-stage base workflow's exact Dura-stripped T1, the
selected Java `AlgorithmMultipleToads` implementation and PyMaCRUISE both
stopped naturally at iteration 43. Their terminal classification arrays are
identical: **0 of 6,092,928 voxels differ**, and Dice is **1.0 for all eleven
labels**. This is selected-source equivalence for the second TOADS stage on
this synthetic workflow input. It does not establish full MATLAB/JIST
composite equivalence or real-MRI validity.

This input differs from the earlier held-fixed replay in the
[retained-center source-fix report](toads-retained-registration-center-source-fix-2026-10-09.md),
where the second pass stopped at iteration 24. The current integrated input
comes from the corrected first TOADS and Dura stages and stops at 43 in both
implementations. The distinct inputs and outcomes should not be conflated.

| Evidence | Value |
|---|---|
| Implementation commit | `b785bff` (`Match retained TOADS registration center`) |
| Current integrated T1 | `build/base-workflow-2026-10-09-retained-center/stages/dura/attempt-01/stripped_t1.nii.gz` |
| T1 SHA-256 | `66a24de1ab12163f786cf6a7e24c2d5707074c3367673f0ced6ec30bca13f0d9` |
| Python classification SHA-256 | `cd0c2a1f058221772b65e9f14d5d6941595e61b0743946c5072913a1dbe7bb3d` |
| Java raw classification SHA-256 | `12ea971a3a4af671754c6d926a9fcdfbd6a2165e351b7546d3ce7d206d5dbdcd` |
| Changed voxels | `0 / 6,092,928` |
| Per-label Dice | `1.0` for labels 1, 5, 6, 14, 15, 16, 18, 22, 24, and 25 |
| Java peak sampled RSS | `3.758 GB` under a `16 GB` watchdog (`-Xmx12g`) |

## Reproduction

The selected Java harness and atlas descriptor are committed/reference build
artifacts in this workspace. The T1 was serialized as little-endian float32,
x-fastest (Fortran order), to match the harness's source-array convention.

```sh
./.venv/bin/python - <<'PY'
import nibabel as nib
import numpy as np
from pathlib import Path

src = Path("build/base-workflow-2026-10-09-retained-center/stages/dura/attempt-01/stripped_t1.nii.gz")
dst = Path("build/base-workflow-2026-10-09-retained-center/java-second/input.raw")
volume = np.asarray(nib.load(str(src)).dataobj, dtype="<f4")
dst.write_bytes(volume.tobytes(order="F"))
PY

CP="$(cat build/java-reference/toads-selected-wrapper-harness/classpath.txt):build/java-reference/toads-selected-wrapper-harness"
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  /usr/bin/java -Xms512m -Xmx12g -cp "$CP" \
  ToadsWorkflowSourceHarness \
  build/base-workflow-2026-10-09-retained-center/java-second/input.raw \
  build/java-reference/source-vandy-atlas.txt \
  172 216 164 0.9375 0.9375 0.9375 50 \
  build/base-workflow-2026-10-09-retained-center/java-second/classification-java.raw

./.venv/bin/python scripts/compare_toads_java_classification.py \
  build/base-workflow-2026-10-09-retained-center/java-second/classification-java.raw \
  build/base-workflow-2026-10-09-retained-center/stages/toads-second/attempt-01/classification.nii.gz
```

The Java run exited successfully, printed `completed=true`, and reached its
natural terminal state at iteration 43. The comparator emitted the per-label
Dice table and zero-difference count above. This is Java/Python stage evidence;
the seven-stage run itself remains synthetic integration evidence.
