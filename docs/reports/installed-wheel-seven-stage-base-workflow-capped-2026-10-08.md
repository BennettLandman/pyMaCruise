# Installed-wheel seven-stage synthetic base workflow (2026-10-08)

## Result and scope

Built `pymacruise-0.1.0-py3-none-any.whl` from implementation commit `a345e2b847db410a5ba2dea33bb2e8f173d12831` and installed it to an isolated target directory with `pip --no-deps`. The command ran from a build directory outside the package checkout with `PYTHONNOUSERSITE=1` and `PYTHONPATH` pointing only to the installed target. The atlas descriptor and atlas maps resolved from that installed wheel.

The installed package completed all seven base stages: N4, first TOADS, Dura, second TOADS, CRUISE core, MaSeg refinement, and surface/thickness. The manifest reports `status=complete`, all seven stages complete, 46 stage-output records, and 24 final-output records. Every stage check passed. The wheel completed the scientific products including the three RAS surfaces, thickness-enriched inner surface, and `.nii.gz` level-set thickness output. The workflow is explicitly marked `converged=false`: both TOADS passes were intentionally capped at two iterations for this packaging check.

An independent inspection verified all 46 output SHA-256 values, all 31 stage NIfTIs as finite and RAS, and all four VTK meshes for finite vertices, triangular faces, and valid point indices. The main run's peak sampled RSS was **6.874 GB** under a 16 GB watchdog. Same-settings resume from the installed wheel hash-verified and skipped all seven stages; peak sampled RSS was **0.123 GB** under a 4 GB watchdog.

This is installed-wheel packaging, bundled-atlas, seven-stage execution, output, stage-check, and restart evidence. The two-iteration TOADS caps mean it does not establish TOADS convergence. It is synthetic evidence and does not establish source-composite equivalence or real-MRI behavior.

Wheel SHA-256: `2821431a2b91e70017eedf149d601ecbd18eca0964bccf79c5afe1c98ee554e5`.

## Reproduction

From the repository root, build and install the wheel:

```bash
mkdir -p build/wheel-base-workflow-current/dist build/wheel-base-workflow-current/site build/wheel-base-workflow-current/tmp
TMPDIR="$PWD/build/wheel-base-workflow-current/tmp" ./.venv/bin/python -m build --wheel --no-isolation --outdir build/wheel-base-workflow-current/dist .
./.venv/bin/python -m pip install --no-deps --target build/wheel-base-workflow-current/site build/wheel-base-workflow-current/dist/pymacruise-0.1.0-py3-none-any.whl
cd build/wheel-base-workflow-current
```

Run from outside the package checkout, pointing Python only at the installed wheel target. The synthetic T1, brain mask, and MaSeg are committed under `tests/data/synthetic/base_workflow_brain_resolution/`:

```bash
PYTHONNOUSERSITE=1 PYTHONPATH="$PWD/site" MPLCONFIGDIR="$PWD/mpl-cache" \
  /Users/landmaba/dev/mipav/pymacruise/.venv/bin/python \
  /Users/landmaba/dev/mipav/pymacruise/scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  /Users/landmaba/dev/mipav/pymacruise/.venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --atlas-descriptor "$PWD/site/pymacruise/data/atlases/cruise-atlas-10obj-toads2008.txt" \
  --maseg /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --output-dir run-current-capped \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 2 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

Append `--resume` to the same `base-run` command and use a 4 GB watchdog to verify that the installed package skips all seven hash-verified stages.

## Evidence and limits

- Installed import was resolved from `build/wheel-base-workflow-current/site/pymacruise/`; no source checkout was on `PYTHONPATH`.
- The manifest's final output products include CRUISE cropped maps and vector field, native-grid level sets and refined labels, three surfaces, thickness outputs, summaries, and provenance.
- Both TOADS stage logs show `iterations=2`, `iteration_limit_reached=true`, and `energy_threshold_reached=false`.
- The earlier [source-checkout brain-resolution report](base-workflow-seven-stage-brain-resolution-2026-10-08.md) provides natural TOADS convergence evidence (iterations 12 and 34) for this same synthetic workflow. It is separate from this capped installed-wheel run.
- Selected Java component comparisons and known deviations are recorded in their linked component reports; neither workflow report claims whole-workflow MATLAB/JIST equivalence.
