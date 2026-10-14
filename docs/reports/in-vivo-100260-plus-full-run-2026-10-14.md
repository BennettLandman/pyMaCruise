# Human PLUS integration run: case 100260

**Evidence category:** real-human workflow execution and quality audit. This is not a Java/MATLAB source-equivalence result; a subject-specific PLUS reference output was not supplied.

## Inputs and implementation

The run used the case 100260 T1-weighted image, high-resolution FLAIR, matching T1 brain mask, and SLANT MaSeg labels. It used the bundled 12-object LesionTOADS descriptor and matching topology atlas. Input data stayed outside Git. The implementation was at commit `ffc40355` (`Allow signed interpolation values in N4 adapter`). The signed-value N4 change was required because registered FLAIR interpolation contains negative values; the ANTsPy N4 adapter completed on that data and returned finite output.

## Reproduction command

From the `pymacruise` repository root, with the project environment and input paths available:

```bash
MPLCONFIGDIR=/private/tmp/mpl-pymacruise PYTHONNOUSERSITE=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src \
  ./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-mp2rage_T1w.nii.gz \
  --flair /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-highres_FLAIR.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-mp2rage_T1w_bet_mask.nii.gz \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/100260/sub-010260/ses-01/SLANT-TICVv1.2/post/FinalResult/sub-010260_ses-01_T1w_seg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/in-vivo-100260-plus-lesionatlas-2026-10-14 \
  --regrid-inputs
```

The first attempt used the BASE 10-object descriptor and failed at LesionTOADS initialization because its class table did not match the PLUS input. The corrected run used the bundled 12-object LesionTOADS descriptor above. It completed all nine stages: preprocessing, crop, initial alignment, registered recrop, lesion workflow, Dura outputs, native Dura transfer, MaSeg cortical separation, and PLUS surface workflow. LesionTOADS stopped by natural convergence at iteration 13 with energy ratio `8.090717165032402e-05`, below the configured `1e-4` criterion. The complete run's peak sampled RSS was 9.924 GB under a 16 GB watchdog.

The identical command with `--resume` verified the completed run and exited successfully; peak sampled RSS was 1.813 GB under an 8 GB watchdog. An integrity scan checked 93 manifest file references with no missing files or hash mismatches. All 67 NIfTI outputs were finite and RAS.

## Ribbon quality result

The audit used `scripts/audit_plus_ribbon_quality.py --maseg-profile macruise`, which interprets MaSeg values through the source-derived GM and WM label groups in `pymacruise.cortical_separation`. The configured minimum tissue support is 97%; the ribbon target is 5 mm with a ±1 mm median/mean tolerance and bounded surface percentiles.

| Check | Result | Criterion |
|---|---:|---:|
| WM inside inner surface | 97.63% | ≥97% — pass |
| GM between inner and outer surfaces | 43.57% | ≥97% — fail |
| Inner surface outside central level set | 41 voxels | 0 — fail |
| Surface connected components | 127 per surface | 1 per surface — fail |
| Mesh topology | Watertight; Euler characteristic 86 per surface | Connected genus-zero surface — fail |
| Closest-point surface thickness | median 1.064 mm; p05 0.243 mm; p95 2.520 mm | median 4–6 mm; p05 ≥2 mm; p95 ≤8 mm — fail |
| Level-set thickness | mean 2.204 mm | 4–6 mm — fail |

The audit report is `build/in-vivo-100260-plus-lesionatlas-2026-10-14/quality-audit-macruise-profile.json` (human data and generated products remain in ignored `build/`). The audit correctly reports `passed: false`. Thus this run demonstrates restartable PLUS execution and natural lesion convergence, but its current surfaces do not satisfy the declared ribbon-quality requirements. The next integration work should trace the GM/support and disconnected-component failures through MaSeg cortical separation, registration/cropping, and TGDM input fields before making algorithm changes.

## Axis-only RAS comparison

Because the historical human examples use axial axis reorientation in the BASE path, the same PLUS inputs were also run with `--no-regrid-inputs`. This is a controlled geometry comparison, not evidence that the historical PLUS example used this setting. The warning records that FLAIR obliquity is discarded while voxel axes are permuted/flipped without interpolation. The workflow implementation remained the same as `ffc40355`; the label-aware audit and status documentation were updated in `c7199912`.

The second run completed and same-settings resume verified all nine stages. LesionTOADS converged naturally at iteration 13 with ratio `5.5984768550843e-05`; peak sampled RSS was 10.041 GB under a 16 GB watchdog. An integrity scan found 88 manifest references with no missing files or hash mismatches, and all 63 NIfTI outputs were finite and RAS. Resume peak sampled RSS was 2.247 GB under an 8 GB watchdog.

| Diagnostic | Default regrid | Axis-only `--no-regrid-inputs` |
|---|---:|---:|
| LesionTOADS terminal iteration / ratio | 13 / `8.090717e-05` | 13 / `5.598477e-05` |
| WM-inside support | 97.63% | 97.63% |
| GM-between support | 43.57% | 43.49% |
| Components per inner/central/outer mesh | 127 / 127 / 127 | 121 / 121 / 121 |
| Inner-surface nesting violations | 41 voxels | 46 voxels |
| Surface median thickness | 1.064 mm | 1.063 mm |
| Level-set mean thickness | 2.204 mm | 2.207 mm |

The geometry choice changes some details but does not explain the very low GM ribbon support or thickness. Both modes complete and converge; both anatomy audits flag the same quality concerns. Without a subject-specific PLUS source run and settings, these are diagnostic observations for locating differences rather than source-equivalence failures.

## Focused validation

After the signed interpolation fix, `tests/test_n4.py` passed **3 tests**. The full pytest result in the README predates commit `ffc40355`; no new full-suite result is claimed here. The MaCRUISE-label audit completed under an 8 GB RSS watchdog with 1.011 GB peak sampled RSS. A synthetic-profile audit remains separate evidence and does not substitute for this human case or for source equivalence.
