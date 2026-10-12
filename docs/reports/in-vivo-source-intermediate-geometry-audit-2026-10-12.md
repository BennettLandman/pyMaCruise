# Human source intermediate geometry audit — 2026-10-12

## Scope

This is a diagnostic comparison of the downloaded Singularity intermediate files with the completed PyMaCRUISE no-regrid run. It does not establish source equivalence. Source stage NIfTIs are not consistently stored on the raw T1 grid: the N4 output has the oblique RAS affine, while TOADS and Dura outputs use LPS axis codes and some omit the T1 translation/obliquity. To avoid a direct LPS/RAS array-order error, the comparison below first converts both NIfTIs to closest-canonical RAS, then compares arrays by index. Because the source headers do not preserve a common physical grid, these values are diagnostic only; they cannot by themselves establish anatomical alignment.

## Results

All files have shape `184×256×256`. The source N4 image and Python N4 image both report oblique RAS and matching affines. Their intensity correlation is **0.996625**, with MAE **1,153.48** and maximum absolute difference **36,488.51**. This is strong intensity-pattern agreement but not numerical equivalence.

The TOADS classification comparisons use the orientation-normalized array indices. Each has 11 labels including background; ten non-background labels fall below Dice 0.97.

| Intermediate | Source orientation | Python orientation | Changed voxels after RAS axis normalization | Lowest non-background label Dice |
|---|---|---|---:|---:|
| TOADS first classification | LPS | RAS | 509,212 | 0.216619 |
| TOADS second classification | LPS | RAS | 285,212 | 0.414184 |

The second-pass classification is closer than the first by changed-voxel count, but both differ substantially. Selected tissue-map diagnostics are also mixed: TOADS-first WM-fill correlation is 0.977736, GM correlation 0.936271, and WMMask Dice 0.928565; TOADS-second WM-fill correlation is 0.987652, GM correlation 0.939576, and WMMask Dice 0.951852. These scalar-map correlations do not replace label-specific acceptance.

To check the axis normalization against source behavior, the selected JIST `ReorientVolume.solve(..., AXIAL, UNCHANGED, WINDOWED_SINC)` was replayed on the downloaded source N4 NIfTI with the supplied MIPAV/JIST runtime. The output remained `184×256×256`, at `1 mm` resolution, with MIPAV axis orientation `[R2L, A2P, I2S]`. Its voxel array compared with the source N4 array after x/y orientation flips at correlation **0.999998**, MAE **125.36**, and maximum absolute error **2,159.25** intensity units. This supports the RAS/LPS array-order normalization used above. It does not validate SPECTRE's atlas mask, stripped image, affine atlas registration, or the mapping of later outputs back to the oblique raw-T1 world coordinates. The replay is reproducible with:

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/compare_human_spectre_reorient.py \
  --output-dir build/invivo-base-oblique-no-regrid-2026-10-12/spectre-reorient-replay
```

The Dura-stripped T1 image has correlation **0.978367**, MAE **924.20**, and maximum absolute difference **194,617.75** after axis normalization. The source and Python Dura brain masks have Dice **0.935051**, with 177,375 differing voxels. Their stored headers also differ in orientation and physical-grid representation.

## Interpretation and next step

The new evidence locates substantial discrepancies before CRUISE/TGDM and MaSeg refinement, so the final surface mismatch cannot be assigned to the CRUISE core alone. The supplied VandyCRUISE source clarifies the intended preprocessing: `ToadsCruiseLong.java` sets SPECTRE “Output in Original Image Space” and “Resample data to Isotropic voxels” to false (lines 91–92). `MedicAlgorithmSPECTRE2010.java` then calls `ReorientVolume.solve(..., AXIAL, UNCHANGED, WINDOWED_SINC)` and sets the output axis directions to R2L/A2P/I2S (lines 376–397); its later `outputInOriginalSpace()` call is commented out (line 695). Thus historical TOADS does run in SPECTRE's axial, LPS-oriented space. The selected axial reorientation replay supports x/y orientation flipping of the downloaded N4 data, but the downloaded run still lacks the compiled wrapper settings and reliable spatial headers for later SPECTRE/TOADS intermediates. The next comparison should recover the later transfer/header mapping well enough to report physical-space source metrics. This is diagnostic only: the selected Python runtime intentionally stays in native T1 RAS and omits BET/SPECTRE, as specified in `docs/toads-port-plan.md`; do not add SPECTRE to that runtime as part of this mapping task. Keep existing index-normalized scores as diagnostic baselines, not physical-space agreement.

Source files reviewed: `toads-cruise/src/edu/vanderbilt/masi/plugins/CRUISE/ToadsCruiseLong/ToadsCruiseLong.java`, `toads-cruise/src/edu/jhu/ece/iacl/plugins/segmentation/skull_strip/MedicAlgorithmSPECTRE2010.java`, and `toads-cruise/src/edu/jhu/ece/iacl/algorithms/volume/ReorientVolume.java` in the supplied `toads-cruise-scm-2026-09-26` tree.

The comparison was generated with NiBabel `as_closest_canonical` on each source/Python pair, same-shape checks, and per-label Dice or scalar MAE/correlation. Source paths are under `sample-macruise/output/ToadsSeperate`; PyMaCRUISE outputs are under the ignored `build/invivo-base-oblique-no-regrid-2026-10-12/stages/` tree. No human data were added to this repository.
