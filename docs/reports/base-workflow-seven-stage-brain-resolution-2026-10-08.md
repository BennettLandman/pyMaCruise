# Seven-stage synthetic base MaCRUISE workflow — 2026-10-08

## Scope

This report records a full brain-resolution synthetic run of the restartable base workflow after integrating its scientific surface and thickness products. It validates that the seven implemented stages compose, converge where configured, write the expected RAS products, pass stage checks, and resume from recorded hashes. It is synthetic integration evidence; it does not claim whole-workflow equivalence to MATLAB/JIST or validate a human MRI.

## Reproduction

Implementation source fingerprint: `b200db7a08b9bec962249dab5de65aa0be11f650324a8e8699b65e42256a748f`. The source was committed in `1b01838`; the full suite had already passed on the same code at `5a412b7` (subsequent commits were documentation only).

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/ace-indexed-heap-2026-10-08/base-seven-stage-run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

Exit code was 0; all seven stage checks passed. The watchdog sampled a peak of **8.433 GB RSS** and did not terminate the run. Both TOADS passes reached their energy thresholds rather than their iteration caps:

| Stage | Iterations | Terminal energy ratio | Stop |
|---|---:|---:|---|
| First TOADS | 12 | `6.96809075e-05` | Energy threshold reached |
| Second TOADS | 34 | `5.54052858e-05` | Energy threshold reached |

The TOADS iteration counts and terminal ratios match the previously recorded selected-Java comparisons for these inputs. This is selected component evidence, not a complete source-composite comparison.

## Outputs and stage checks

The manifest at `build/ace-indexed-heap-2026-10-08/base-seven-stage-run/run.json` reports `status=complete`, `converged=true`, and the scope `synthetic/native-grid base route; not a full Java-equivalence claim`. Its seven complete stages are N4, first TOADS, Dura, second TOADS, CRUISE core, MaSeg refinement, and surface/thickness. All **46 stage output records** are hash-checked. The 24 final-output records include cropped CRUISE tissue/ACE/GVF products, full-grid level sets, refined MaSeg labels and summary, three RAS-mm VTK surfaces, an inner surface with closest-point thickness fields, a finite level-set thickness `.nii.gz`, summaries, provenance, and logs.

An independent post-run inspection loaded all **31 stage NIfTIs** and confirmed RAS orientation and finite voxel values. It parsed all four VTK files and confirmed finite RAS-mm points, valid triangular faces, and in-range point indices. The core CRUISE crop and vector field are on `(166,208,156)`; the inner, central, outer, refinement, and thickness volumes are on the native `(172,216,164)` grid. The thickness summary reports 142,858 crossing voxels, level-set mean **2.00259464 mm**, and symmetric closest-point surface mean **1.64906598 mm**. Surface extraction remains the documented scikit-image Lewiner substitution for JIST `IsoSurfaceOnGrid` connectivity 18.

## Restart check

The same command and settings, with `--resume`, exited 0 and verified/skipped all seven stages. Under a 4 GB watchdog, peak sampled RSS was **0.124 GB**. No stage was recomputed.

## Test evidence and limits

- Focused restart/integration test: `./.venv/bin/python -m pytest -q tests/test_base_workflow_restartable.py` — **3 passed in 32.43s** on implementation commit `1b01838`.
- Full suite: **1,243 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 1,743.86s**, peak sampled RSS **14.961 GB**, on `5a412b7`; see [full-suite report](full-pytest-2026-10-08-base-surface-thickness.md). This full suite ran before the full-grid integration, but on the same implementation code.
- This run does not establish end-to-end MATLAB/JIST equivalence. Remaining base-workflow source comparisons include MaSeg refinement and composite outputs; Lewiner/JIST mesh topology remains a known substitution risk. Plus-workflow parity and real-MRI validation also remain open. The available subject T1 still lacks its matching user-provided brain mask.
