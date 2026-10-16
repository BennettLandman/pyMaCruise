# Exact-input sample MaCRUISE BASE run — 2026-10-16

## Result and evidence category

PyMaCRUISE completed and resumed the full eight-stage BASE workflow on the
184×256×256 T1 and MaSeg that accompany the downloaded Singularity output. The
source `Raw/T1.nii.gz` and `Raw/Seg.nii.gz` voxel arrays are byte-identical to
the supplied input arrays (SHA-256 `aae1e316…068b` and `524e0881…e81`). The
run supplied the matching T1-grid BET brain mask, used the V3 BASE settings
listed below, and followed the requested axis-only RAS mode: voxel values are
not interpolated and the oblique affine directions are discarded with
warnings.

This is a **real-human workflow execution plus a diagnostic comparison to a
historical output tree**. It is not controlled whole-workflow equivalence:
the Singularity output directory has no run log or manifest, so its exact
source revision, selected branch, settings, and upstream inputs beyond the
raw T1/MaSeg arrays are unverified. The historical output includes results of
its own skull-stripping path; PyMaCRUISE follows the owner-selected
user-provided-mask contract.

## Reproduction

Run from the repository root:

```bash
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/sample-singularity-base-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/input/T1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --output-dir build/sample-singularity-base-source-profile-2026-10-15 \
  --no-regrid-inputs --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --dura-membership-level 0.3 \
  --crop-edge-width 5 --ace-diffusion-iterations 25 --gvf-iterations 200 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 --wm-iso-level 0.5 \
  --ace-use-updated-csf-in-ace 2 --source-faithful-refinement
```

Resume with the same options and append `--resume`. It verifies the run
signature and every completed stage before skipping work.

## Execution and integrity

- All **8/8 workflow stages** completed. Both TOADS passes stopped naturally:
  pass 1 at iteration 20, ratio `4.29838547e-05`; pass 2 at iteration 18,
  ratio `4.44517791e-05`; the configured threshold was `1e-4` and cap 50.
- Peak sampled RSS was **5.420 GB** under a **16 GB** watchdog; it did not
  trigger.
- Resume verified all eight stages and exited successfully, peaking at
  **0.699 GB** under an 8 GB watchdog.
- Independent manifest verification found **68/68 output hashes valid** and
  **51/51 NIfTIs finite and RAS**. The output run manifest and products are
  retained under ignored `build/sample-singularity-base-source-profile-2026-10-15/`.

## Controlled selected-Java TOADS handoffs

To separate TOADS implementation parity from the uncertain historical run, the
selected VandyCRUISE Java wrapper was run on each exact intensity array
prepared for Python TOADS. For each pass, its input is the same float32 array
with voxels outside the supplied brain mask set to zero; dimensions are
184×256×256, spacing is 1 mm, atlas descriptor is the bundled
`cruise-atlas-10obj-toads2008.txt`, and the selected-wrapper iteration ceiling
is 50. Java's raw result is reshaped in source-compatible Fortran order and
compared voxel-for-voxel with the Python stage classification.

- First handoff: **0/12,058,624 changed voxels**; all 11 labels have Dice
  **1.0**. Java completed successfully under a 16 GB RSS watchdog with peak
  sampled RSS **3.299 GB**.
- Second handoff: **0/12,058,624 changed voxels**; all 11 labels have Dice
  **1.0**. Java completed successfully under the same watchdog with peak
  sampled RSS **4.409 GB**.

These are controlled, held-fixed **TOADS-stage source equivalence results on
this human-derived data**, not whole-workflow equivalence: Python generated
the N4/Dura handoff arrays and the complete original workflow provenance is
still unavailable. Machine-readable per-label results are in the run's
ignored `build/sample-singularity-base-source-profile-2026-10-15/java-first/comparison.json`
and `.../java-second/comparison.json`; the raw Java inputs and outputs are
also retained in ignored `build/` for local reproduction. The source harness
is `tests/java_reference/ToadsWorkflowSourceHarness.java`; invoke it with the
bundled source classpath, the 10-object atlas descriptor, dimensions,
spacing, 50 iterations, and each handoff's raw input/output paths.

## Historical-output comparison

The comparison first verified exact source/Python input voxel identity, then
reoriented source outputs by axis permutations/flips only and compared arrays
by voxel index. Python mesh vertices were mapped through the Python voxel
grid and the source `Raw/T1` affine before measuring distances in millimeters.
This honors the no-interpolation comparison geometry but cannot recover the
historical workflow settings.

- Cached source-vs-input brain-mask Dice is **0.95797**; cached Dura brain-mask
  Dice is **0.96790**. The source WM and atlas-mask products match Python
  exactly; source GM and brain-tissue masks each differ by 13,885 voxels
  (Dice **0.98908** and **0.99354**).
- First TOADS classification differs at **509,194/12,058,624** voxels; mean
  nonzero-label Dice is **0.80694** (minimum **0.21580**). Second-pass
  classification differs at **198,646** voxels; mean label Dice is
  **0.88173** (minimum **0.59263**). These are comparisons to saved outputs,
  not same-input selected-Java comparisons.
- N4 corrected intensity MAE is **6,165.91** on the same indexed array. The
  Dura-stripped T1 and Dura mask also differ substantially (T1 MAE
  **3,859.87**, mask Dice **0.96790**); therefore differences are already
  present upstream of the cortical surface stages.
- Historical refinement has **37,372** differing voxels across 134
  nonbackground labels; mean label Dice is **0.97262**, minimum **0.44**, and
  **35 labels** are below the agreed 0.97 minimum. This comparison does not
  override the fresh-MATLAB result: fresh MATLAB matches the corrected Python
  refinement exactly, while the provenance of this saved historical result
  remains unknown.
- Inner, central, and outer surface mean absolute distances are **0.1441,
  0.2189, and 0.2664 mm**; sampled Hausdorff distances are **4.138, 5.831,
  and 7.059 mm**. None meets the agreed ≤0.01 mm mean / ≤1 mm Hausdorff
  limits. These mismatches cannot be attributed to surface extraction alone;
  upstream maps differ too.

The Python anatomical diagnostic reports GM-between-surface support **85.23%**
and WM-inside-inner support **88.42%** on the supplied MaSeg, with 4 inner
voxels outside central, 0 central voxels outside outer, and 132 inner surface
voxels outside the supplied mask. Meshes are watertight but have 67
components each (Euler characteristic −198); level-set mean thickness is
4.226 mm and symmetric closest-point surface thickness mean is 2.653 mm.
These measurements describe this run and are not substitutes for source
equivalence gates.

## Interpretation and next integration work

This run establishes a restartable full BASE execution on the exact raw image
arrays present in the downloaded output and exposes the first differences at
N4/mask selection, then TOADS and Dura. The failed historical comparisons do
not establish that the Python implementation is incorrect because the source
run provenance is missing and its skull-stripping/input-mask contract differs.
For a controlled decision, capture a fresh original V3 BASE run on these
voxel arrays with the supplied mask replacing skull stripping, recording the
compiled VandyCRUISE/JIST versions, all wrapper settings, intermediate
products, and hashes. Compare those outputs stage by stage before changing
the Python algorithms.
