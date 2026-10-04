# Brain-resolution source-order TOADS through 12 iterations

## Run and result

The production `pymacruise toads-loop` command completed on the deterministic
synthetic `(172,216,164)` T1/brain-mask fixture and the bundled Vandy atlas. The
run reached the requested 12-iteration cap but did **not** reach the `1e-4`
energy threshold: terminal ratio was **`0.0132785104`**. The last `stop=true`
means iteration limit, not convergence (`energy_threshold_reached=false`).

Command from the repository root:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise toads-loop \
  --t1-n4 build/vandy-brain-resolution-memory/t1.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output build/vandy-brain-resolution-memory/toads-source-order-12-iterations.nii.gz \
  --maximum-iterations 12 --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --padding 10 --relation-limit 3 --maximum-distance 1e-4 \
  --provenance build/vandy-brain-resolution-memory/toads-source-order-12-iterations.json
```

The command exited 0 in **1,845.5 seconds**. Peak sampled RSS was **7.995 GB**
under the 12 GB watchdog; it was not terminated. The classification is finite
uint8 RAS NIfTI on the native `(172,216,164)` grid with labels
`[1,5,6,14,15,16,17,18,22,24,25]`. SHA-256:
`d9a40ac33706b68e8f5a136f1133933ff54183a52af7c28246842025cd6720e4`.
Outputs and full run manifest are ignored build artifacts.

| Iteration | Energy ratio | Stop reason |
|---:|---:|---|
| 1 | 1.66407061 | continue |
| 2 | 0.323834717 | continue |
| 3 | 0.00333001465 | continue |
| 4 | 0.0130230254 | continue |
| 5 | 0.485990405 | continue |
| 6 | 0.278568447 | continue |
| 7 | 0.162256882 | continue |
| 8 | 0.0381764956 | continue |
| 9 | 0.0127420947 | continue |
| 10 | 0.00365858758 | continue |
| 11 | 0.0164114311 | continue |
| 12 | 0.0132785104 | maximum iteration cap |

## Interpretation

The first two ratios match the prior source-order brain-resolution TOADS run,
confirming the synthetic input and startup path. Later ratios are
non-monotonic: a sharp drop at iteration 3 is followed by a large increase at
iteration 5, then a general reduction with further rebounds. The energy
threshold remains unreached. This is meaningful Python workflow evidence, not
Java equivalence; the compact Java/Python comparison still diverges after
iteration 2. The supplied human MRI was not used.

Next, localize the first later-state difference against Java on a compact
source-order fixture, then determine whether the same corrected state updates
improve brain-resolution convergence before running the complete base
workflow again.
