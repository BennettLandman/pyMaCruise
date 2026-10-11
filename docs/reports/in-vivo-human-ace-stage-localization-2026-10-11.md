# Human MaACE stage localization — 2026-10-11

## Scope

This is a held-fixed, real-human-data source localization using the already
saved corrected boundary-enabled base-tail inputs. It did not rerun tissue
preparation, TOADS, the base workflow, or Java. The common RAS grid is
141×173×131 at 1 mm. The Java intermediate products were captured during the
earlier direct Java MaACE replay, with a 4 GB Java heap and 16 GB RSS watchdog.
Python stage replays also ran under the 16 GB watchdog; their highest sampled
RSS was 0.945 GB. Input and source hashes and wrapper settings are in the
[held-fixed comparison JSON](../../build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/comparison.json).

No general scalar-field ACE tolerance has been declared. The continuous-field
numbers below localize where the paths first differ; they are not a claim that
every reported field passes a pre-agreed numerical tolerance.

The stage table below is the **pre-correction diagnosis**. A source-ordered
FMM correction and exact fixed-input human replay are recorded in the
[follow-up report](in-vivo-human-ace-ordered-fmm-correction-2026-10-11.md).

## Earliest divergence and stage trace

| Java intermediate | Comparison | Result |
|---|---|---|
| `01_tmpVolfirst` | Initial `126.9 - WM` float field | Exact, 0/3,195,483 voxels differ |
| `02_PHIinit` | Initial unit-speed signed distance | First divergence: 602,421 voxels differ; MAE 2.1245e-5, max 0.0282317. The sign mask remains the same. This remains different with Python's optional per-update float32 rounding enabled. |
| `03_csf_naive` | Integer fuzzy-CSF initialization | Exact |
| `04_wm_morph`, `05_closedGM_morph` | Combined WM threshold and four 18-neighbor closing iterations | Both exact |
| `07_tmpVol_diffusion` | Gaussianized signed-distance guide | 2,503,191 values differ; MAE 2.1849e-5, max 0.0074396 |
| `08_initcsf_beforediffusion` | Initial CSF seed field | Exact |
| `09_initcsf_diffusion` | 25 source-faithful diffusion iterations | 1,793,016 values differ; MAE 2.7363e-7, max 0.00072664 |
| `13_PHI_narrowBand` | Narrow-band evolution and reinitialization | 534,886 values differ; MAE 1.1537e-5, max 0.25510; 2 sign-mask voxels differ |
| `14_tmpVol_narrowBand` | CSF-weighted travel distance | 2,492,055 values differ; MAE 1.3140e-5, max 0.25795 |
| `15_csf_beforeclean` | Thresholded ACE skeleton memberships | 206 membership intensities differ, but only 3 foreground voxels differ |
| `16_csf_afterclean` | Component cleanup | Same 3-voxel foreground difference remains |

Thus the first numerical divergence is the initial distance field, while the
first changed ACE seed support appears much later, at the magnitude threshold
after narrow-band evolution and weighted-distance propagation. The initial
CSF seed support, WM morphology, and cleanup are not the cause.

## The three pre-thinning support changes

The source uses `magnitude < 0.85` and then truncates `254 * (magnitude +
0.05)` into the candidate membership. At the two enhanced-GM disagreements:

| Crop voxel | Input GM | Java magnitude → membership | Python magnitude → membership | Effect |
|---|---:|---:|---:|---|
| `(51,34,71)` | 247.4375 | 0.849898 → 228 | 0.850020 → 0 | Java removes this GM voxel; Python retains it. |
| `(110,115,90)` | 236.79486 | 0.856149 → 0 | 0.792624 → 214 | Java retains this GM voxel; Python removes it. |

The third support change is `(52,158,80)`: GM is zero; Java magnitude is
0.850018 (excluded) and Python magnitude is 0.849965 (membership 228). It
changes the skeleton but not enhanced GM. Java and Python speed values at
these points are very close; the candidate flips arise at the magnitude
threshold after their distance/gradient fields diverge.

This resolves why the enhanced-GM aggregate MAE is tiny despite two large
single-voxel errors. It also explains the thinned-skeleton disagreement:
after matching Java's cleaned membership support exactly as the input,
Python `thin_ace` returns the exact Java thinned output (94,941 foreground
voxels; Dice 1.0). With the actual Python ACE seed, the output has Dice
0.941188 against Java and differs at 11,165 voxels. This isolates the
thinning mismatch to the three changed input support voxels, rather than a
separate ThinACE parity problem on this case. The source's distance-priority
simple-point peeling is topology-sensitive, so a few support changes can
alter many later deletions.

The cleaned ACE skeleton before thinning has binary Dice 0.9999899, above the
user's 0.97 segmentation criterion. The thinned skeleton remains below that
criterion. These are selected-stage comparisons on one human scan, not
whole-workflow equivalence or clinical validation.

## Conclusion after correction

The initial mismatch came from the Python unit-speed FMM's duplicate-entry
`heapq` traversal and insertion order, which did not preserve Java's indexed
float-priority heap behavior. `_signed_distance_fmm` now reuses the validated
source-ordered unit-speed kernel; on the saved human initial level set it
matches Java `02_PHIinit` exactly, and the held-fixed ACE products now match
Java exactly. See the follow-up report for the focused tests, replay command,
and resource measurements. The selected Java sources are
[`DistanceFieldyh.java`](../../tests/java_reference/upstream/edu/vanderbilt/masi/algorithms/CRUISE/ace/DistanceFieldyh.java),
[`BinaryMinHeap.java`](../../tests/java_reference/upstream/edu/jhu/ece/iacl/jist/structures/data/BinaryMinHeap.java),
and [`MaAnatomicallyConsistentEnhance.java`](../../tests/java_reference/upstream/edu/vanderbilt/masi/algorithms/CRUISE/ace/MaAnatomicallyConsistentEnhance.java).
The Python operations are in
[`ace.py`](../../src/pymacruise/ace.py#L580) (initial FMM),
[`ace.py`](../../src/pymacruise/ace.py#L806) (narrow band), and
[`ace.py`](../../src/pymacruise/ace.py#L889) (weighted FMM); thinning is in
[`skeleton.py`](../../src/pymacruise/skeleton.py#L70).

## Local reproducibility

All scripts and generated stage arrays are under the ignored
`build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/`
directory. The early-stage replay and the narrow-band-to-candidate replay
were run under the 16 GB watchdog. The Java outputs and intermediate arrays
were reused without rerunning the JVM. The thinning cross-check ran Python
`thin_ace` on the saved Java pre-thinning membership. No pytest or full/human
workflow run was performed for this diagnostic, and no production code,
README, or plan file was changed.

The bounded Python replays can be repeated from the repository root:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/stage-localization/trace_early.py

PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/stage-localization/trace_to_candidates.py

PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/stage-localization/thin_crosscheck.py
```
