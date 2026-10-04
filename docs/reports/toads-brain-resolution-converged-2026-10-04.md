# Brain-resolution source-order TOADS convergence

## Result

The production `pymacruise toads-loop` completed the selected standalone
source-order TOADS route on the deterministic `(172,216,164)` synthetic
T1/brain-mask fixture and bundled Vandy atlas. Configured maximum iterations
were 50; the selected energy criterion was reached at **iteration 20**, with
terminal ratio **`6.04866409e-6`** (threshold `1e-4`). The loop reported
`energy_threshold_reached=true` and `iteration_limit_reached=false`.

Command from the repository root:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  ./.venv/bin/pymacruise toads-loop \
  --t1-n4 build/vandy-brain-resolution-memory/t1.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output build/vandy-brain-resolution-memory/toads-source-order-50-iterations.nii.gz \
  --maximum-iterations 50 --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --padding 10 --relation-limit 3 --maximum-distance 1e-4 \
  --provenance build/vandy-brain-resolution-memory/toads-source-order-50-iterations.json
```

The run exited 0 after **3,515.5 seconds**. Peak sampled RSS was **8.056 GB**
under the 12 GB watchdog; it was not terminated. Output SHA-256:
`45dc28319ba63499fdba8349c55656e4b30d3d155b54be83030feaf6a2b6ce9d`.
The result is finite uint8 RAS NIfTI on the `(172,216,164)` input grid with
labels `[1,5,6,14,15,16,17,18,22,24,25]`; provenance records the atlas and
effective parameters. The image and manifest are ignored build artifacts,
regenerated from the committed synthetic fixture generator.

| Iteration | Energy ratio | Iteration | Energy ratio |
|---:|---:|---:|---:|
| 1 | 1.66407061 | 11 | 0.0164114311 |
| 2 | 0.323834717 | 12 | 0.0132785104 |
| 3 | 0.00333001465 | 13 | 0.00214143074 |
| 4 | 0.0130230254 | 14 | 0.00329869217 |
| 5 | 0.485990405 | 15 | 0.00461003929 |
| 6 | 0.278568447 | 16 | 0.00114005851 |
| 7 | 0.162256882 | 17 | 0.0036779223 |
| 8 | 0.0381764956 | 18 | 0.000171895212 |
| 9 | 0.0127420947 | 19 | 0.00207977463 |
| 10 | 0.00365858758 | 20 | 0.00000604866409 |

## Output stability and limits

The converged iteration-20 classification differs from the same Python run
stopped at iteration 12 in 3,536 voxels. Label-wise Dice between those two
Python outputs ranges from **0.982672** (label 18) to above 0.9997 for most
labels; all 11 labels exceed the requested 0.97 threshold. This is an
iteration-stability comparison between two Python outputs, **not** a Java
reference Dice result.

The trajectory is substantially non-monotonic: for example, it rises from
`0.0130230` at iteration 4 to `0.4859904` at iteration 5, then later reaches
`0.0001719` at iteration 18 before rebounding at 19. Java/Python comparison
evidence still diverges after iteration 2 on the compact fixture. Thus this
result establishes synthetic brain-resolution convergence of the Python
standalone TOADS route, not source equivalence, two-pass MaCRUISE completion,
whole-workflow numerical acceptance, or real-MRI validity.

## Next integration step

Run the restartable six-stage synthetic base workflow with the converged
TOADS iteration settings. Confirm that both TOADS passes, Dura, CRUISE core,
and MaSeg refinement complete and pass the manifest/stage checks, then compare
the converged whole-route outputs with Java reference outputs before real-data
validation.
