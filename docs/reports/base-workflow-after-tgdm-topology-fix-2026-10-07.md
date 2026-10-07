# Corrected synthetic base workflow after TGDM topology fix — 2026-10-07

## Reproduction

The end-to-end command ran on the committed `(172,216,164)` RAS synthetic fixture, with a 16 GB RSS watchdog. The exact settings are preserved in `docs/reports/base-workflow-after-tgdm-topology-fix-2026-10-07.json`; this abbreviated command shows the runnable entry point:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/ace-indexed-heap-2026-10-07/base-run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1 --resume
```

The run manifest reports `complete`, `converged=true`, six successful stage checks, 38 hashed stage outputs (30 finite RAS NIfTIs, two JSON files, six logs), and 16 final products. First TOADS converged at iteration 12 (terminal energy ratio `6.96809075e-05`); second TOADS converged at iteration 34 (`5.54052858e-05`). Peak sampled RSS was **8.106 GB**. A second `--resume` invocation verified/skipped all six stages and peaked at **0.127 GB** under a 4 GB watchdog.

## Selected Java comparisons

These are synthetic component/stage comparisons, not a claim of whole-workflow source equivalence.

- First TOADS: selected Java and Python have matching float32 energy ratios through the 12-iteration natural stop, and exact native-grid classification (0 changed voxels).
- Dura: the Java and Python stage received the same held-fixed T1, GM, WM, and WM-mask arrays (input hashes are in the JSON report). The two topology outputs and WM+GM union are exact. After the TGDM connectivity/float32 correction, the final mask differs at 22 boundary voxels; foreground Dice is `0.9999900124`, symmetric mean surface distance `0.000205366 mm`, and Hausdorff `0.9375 mm`. The TGDM zero-level-set signs have zero disagreements. This meets the agreed segmentation and surface tolerances.
- Second TOADS: selected Java was run on the corrected Dura T1 input from this same base run. All **34** energy ratios match at float32 precision through natural convergence; final classification is byte-identical after x-fastest serialization (0 of 6,092,928 voxels changed; Dice 1.0 for each of 11 labels). The Java golden and provenance are committed under `tests/data/synthetic/base_workflow_brain_resolution/`.

Focused validation on implementation commit `9eb84bc`: `tests/test_tgdm.py` **24 passed in 14.85s**; `tests/test_base_workflow_restartable.py` **3 passed in 41.13s**. This cycle did not run the full suite; its latest documented full run predates the implementation commit. The data are synthetic; no real MRI was used. The complete historical V3 workflow, CRUISE/refinement source output comparison, plus workflow source parity, and real-data validation remain open.

Detailed machine-readable metrics, hashes, all 34 Java/Python energy ratios, label counts, and run-manifest values: [JSON evidence](base-workflow-after-tgdm-topology-fix-2026-10-07.json).
