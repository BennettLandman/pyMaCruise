# Human BASE refinement refresh — 2026-10-14

## Scope and result

After correcting Python's interpretation of MATLAB `tabulate` rows, the case
100260 no-regrid BASE run was resumed with the same inputs and source-selected
settings. The workflow rebuilt all eight stages because the implementation
fingerprint changed. All eight completed; the manifest records **68/68 stage
outputs**, and a same-settings `--resume` verified and skipped every stage.
Peak sampled RSS was **6.987 GB** under a 16 GB watchdog. Both TOADS passes
stopped naturally: pass one at iteration 21 (ratio `4.86711651e-05`) and pass
two at iteration 20 (ratio `1.48473373e-05`). The manifest independently
checks **51/51 finite RAS NIfTIs** and **68/68 output hashes**.

The Python source-faithful refinement output was compared against a fresh run
of the unmodified original `MAsegRefineBySurf.m` using the same human MaSeg,
inner, and outer inputs. The arrays are identical at all **12,058,624** voxels;
all 134 label Dice values are 1.0, and none is below 0.97. Python result SHA:
`d3659f46fd46c3d30b8f62970fe07d8b2450913ca330d35848404fe3c2ff5b9f`;
fresh MATLAB SHA: `113677cb60681e7a2fee25ef85e726b80a23089d53de795f7014af7f5c9ee159`.

The historical downloaded full-pipeline `SegRefine/SEG_refine.nii.gz` is a
different artifact. Against fresh MATLAB it differs at **47,306** voxels, has
mean label Dice **0.965159**, minimum Dice **0.382979**, and **56/134** labels
below 0.97. This remains a provenance discrepancy; the historical output is
not substituted for fresh source execution.

The refreshed level sets and tissue preparation reproduce the prior
source-profile run exactly. Independent historical whole-output comparisons
remain outside tolerance: inner/central/outer nonpositive-mask Dice are
0.97240/0.97119/0.97102, while symmetric surface means are
0.15652/0.24067/0.32186 mm and sampled Hausdorff distances are
5.63056/5.83323/7.05927 mm. Refined-MaSeg GM-between-surfaces and
WM-inside-inner support are 89.74% and 96.29%; the corresponding original
MaSeg supports are 84.04% and 88.00%. Each output mesh has 82 connected
components, with 80–81 components smaller than 100 vertices. Thus the direct
refinement source comparison passes; the historical full pipeline still does
not match at the agreed whole-output tolerances. Selected-JIST extractor
evidence on the unchanged fixed level sets remains separate from these
upstream surface-position differences.

## Reproduction

Run from the repository root with the installed project environment. The
command resumes the saved run when present and is restartable; its output
manifest records each stage's settings, hashes, and convergence data.

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 45 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/input/T1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/invivo-base-orientation-only-source-profile-2026-10-12 \
  --no-regrid-inputs --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --ace-diffusion-iterations 25 \
  --ace-use-updated-csf-in-ace 2 --wm-iso-level 0.5 --gvf-iterations 200 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 --resume
```

Detailed, machine-readable artifacts are under the ignored `build/` run
directory: `run.json`, `refreshed-quality-audit.json`,
`refreshed-source-comparison.json`,
`refreshed-refinement-source-comparison.json`, the MATLAB replay/output, and
per-stage logs/products. Human images and original source outputs were not
committed. This is real-human stage/workflow execution plus a direct MATLAB
refinement comparison; it is not whole-BASE source equivalence or a pass of
historical surface-position thresholds.

## Full-suite evidence

The full suite passed on implementation/evidence commit
`f03cdfe2c1299bd929abcd4752a747968be40700`: **1,577 passed, 4 skipped, 0
failed**, 18 expected orientation warnings, 732.76 seconds, and **10.197 GB**
peak sampled RSS under a 16 GB watchdog. Command:

```sh
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 45 -- \
  .venv/bin/pytest -q --tb=short
```

The repository root is included in `PYTHONPATH` because one test imports the
repository's `scripts` namespace. An initial attempt without that path stopped
during collection; it is not the reported full-suite run.
