# Human base tail with the V3 ACE boundary active

## Run and execution evidence

This corrected human tail starts from the completed 2026-10-11 base-run
manifest and reuses its hash-verified `toads-second` products. It runs
base-tissue preparation, ACE/GVF/TGDM, source-faithful MaSeg refinement, and
surface/thickness generation with the MaSeg-derived CSF boundary supplied to
ACE mode 2. The output is
`build/invivo-base-source-v3-human-boundary-2026-10-11/run.json`; inputs and
outputs remain in ignored `build/`.

```sh
.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli base-resume-from-toads \
  --source-run build/invivo-base-source-v3-human-rasgrid-2026-10-11/run.json \
  --maseg build/invivo-base-source-v3-human-inputs-2026-10-11/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --output-dir build/invivo-base-source-v3-human-boundary-2026-10-11
```

The tail completed all four stages: base-tissue preparation (22 output
records), CRUISE core (13), MaSeg refinement (3), and surface/thickness (8).
All stage outputs are hash-recorded; the audit found **33 finite RAS NIfTIs
among 43 files**, with no invalid NIfTI geometry or voxel values. The
12-GiB-watchdog run exited successfully with peak sampled RSS **3.008 GiB**.
A same-settings `--resume` completed in 3.07 seconds with peak RSS **0.752
GiB**; all four stages remained at attempt 1 with their recorded output
hashes unchanged.

The final boundary was created on the shared cropped RAS grid
`(141, 173, 131)`: raw boundary **147,611 voxels**, after second-pass TOADS WM
cleanup **147,282**, and after updated-WM cleanup **140,080**. The CRUISE-core
command and provenance identify this exact `boundary_csf.nii.gz` input, with
ACE mode 2 and 25 diffusion iterations. The selected Java `IsoSurfaceOnGrid`
comparison consumes the resulting fixed level sets and validates extraction
only; see [the comparison report](in-vivo-boundary-jist-surface-extraction-2026-10-11.md).

## Human ribbon and mesh checks

Against the supplied MaSeg, **529,113 of 642,291 GM voxels (82.3790%)** lie
between the inner and outer surfaces, and **379,000 of 428,303 WM voxels
(88.4888%)** lie inside the inner surface. These remain below the requested
0.97 Dice threshold; they are support/containment fractions, not Dice. The
source-faithful refined segmentation keeps the same supported voxel counts
and changes the denominators, producing 89.0175% GM and 96.4747% WM
containment. Three inner voxels lie outside central; none of central lie
outside outer; 304 inner voxels lie outside the supplied brain mask.

All three meshes remain edge-watertight and have 143 connected components,
142 smaller than 100 vertices. Level-set subtraction reports mean thickness
4.1844 mm; mean closest-point surface distance is 2.6253 mm. These describe
this result and are not reference accuracy measurements. The new support
fractions are only 0.0378 and 0.0049 percentage points above the earlier
boundary-free run, respectively; the change is small and does not address the
support failure. No causal conclusion is drawn from this difference.

On the exact corrected level sets, direct selected Java and Python surface
extraction have identical vertex, face, and component counts for inner
(234,234/468,468/143), central (235,178/470,356/143), and outer
(303,098/606,196/143). Symmetric mean vertex-to-triangle distances are
0.00000012/0.00000013/0.00000014 mm; sampled Hausdorff estimates are
0.000345/0.000345/0.000488 mm. These pass the user's 0.01 mm mean and 1 mm
Hausdorff thresholds for the fixed-field extraction component; they do not
validate anatomy, upstream stages, or the entire historical workflow.

Machine-readable local outputs:

- `build/invivo-base-source-v3-human-boundary-2026-10-11/human-quality-audit.json`
- `build/invivo-base-source-v3-human-boundary-2026-10-11/support-localization/`
- `build/invivo-boundary-jist-surface-comparison-2026-10-11/comparison.json`

## Interpretation and limitations

The source boundary-to-ACE wiring is now demonstrated on synthetic phantoms
and exercised in a full human base tail with restart validation. The earlier
human metrics remain historical diagnostics from a run where ACE mode 2 had
no boundary input. Both the earlier and corrected human runs fail the agreed
anatomical support target. The supplied Singularity output directory has no
saved historical stage or final outputs, so full-workflow Java/MATLAB
equivalence cannot be computed from the available data. The run reuses the
verified upstream TOADS checkpoint and does not re-execute registration or
TOADS; this is human integration evidence, not whole-workflow equivalence or
clinical validation.
