# Adam BASE corrections and measured comparison — October 29, 2026

## Execution and engineering checks

The included 0.1.5 wheel was separately installed and passed **128 targeted
tests**, with 14 orientation warnings, in
167.07 seconds. The full project suite was not rerun.
All ten stages of the independent raw-input Python reconstruction completed.
Both TOADS passes reached their energy thresholds naturally at iterations
21 and
24. The supplied checker
verified 84 stage artifacts, all final image grids, integer labels,
finite products and mesh integrity. These checks establish engineering
completion; anatomical comparison is reported separately below.

This case validation uses explicitly audited checkpoint continuation. The
first five stages were independently computed from the supplied raw T1, MaSeg
and brain mask with frozen 0.1.2 Python code. The final five stages used this
0.1.5 wheel. Wheel comparison verifies identical upstream CLI branches,
upstream orchestration and shared downward-topology kernels. Upstream inputs,
assets, derived resources and outputs were hash-verified before continuation;
upstream outputs were rechecked afterward. The added upward topology branch
is selected only by downstream BASE core. Shared default behavior is covered
by the regression tests. This is a complete independent Python chain, not a
monolithic fresh 0.1.5 invocation on this human scan. The installed-wheel
synthetic BASE tests exercise the complete main command and restart behavior.
No original MaCRUISE intermediate was fed into the reconstruction.

All five final NIfTIs are **256×256×150**, exactly match the original T1 affine,
and preserve qform/sform and storage axes. SEG_refine is int32 with zero
fractional voxels. Primary surfaces use physical RAS mm. The optional bundled
comparison tool creates explicitly named voxel-index copies without changing
the primary outputs or fitting a registration.

## Independent reconstruction comparison on Adam's scan

The before and after chains start from the same supplied raw inputs. Original
products are used only for comparison. Original raw VTK coordinates are mapped
through the original Raw/Seg affine; this mapping matches the supplied FreeView
ASC exports to below 0.000009 mm. No fitted registration, smoothing, component
filtering or cosmetic label cleanup is used to improve agreement.

| Label measure | 0.1.1 | 0.1.5 |
|---|---:|---:|
| Differing voxels | 97,837 | 125,524 |
| Nonzero labels below Dice 0.97 | 84 | 61 |
| Mean per-label Dice | 0.93822 | 0.93193 |

Differing label voxels increased **28.3%** relative to 0.1.1.
Central and inner mean distances also worsen, while outer mean distance and
component counts improve. This mixed result is not an overall accuracy
improvement. Preserve the previous release and compare both outputs.
All 134 nonzero label comparisons are retained in
`validation/comparison.json`.

| Surface | Mean mm before | Mean mm after | Sampled max mm before | Sampled max mm after | Components before → after |
|---|---:|---:|---:|---:|---|
| Inner | 0.36073 | 0.55236 | 7.73516 | 10.17665 | 36 → 1 |
| Central | 0.34652 | 0.42910 | 8.27032 | 9.74151 | 36 → 1 |
| Outer | 0.26585 | 0.17433 | 8.22603 | 11.20379 | 36 → 1 |

The mean averages both directed all-vertex-to-opposite-triangle means. The
sampled maximum is not a continuous triangle-to-triangle Hausdorff bound.
Original meshes have one component each. Project agreement gates require
each label Dice ≥0.97, symmetric mean surface distance ≤0.01 mm and sampled
maximum ≤1 mm. Combined result: **False**.
Remaining differences must not be described as fully resolved.

## Source-derived fixes and isolated validation

1. **Topology propagation:** the supplied v3.2 `run_TopoCorrect_jist.m` selects
   intensity initialization, background-to-object propagation and disabled
   intensity threshold stopping. Python previously used paint-mask propagation
   in the opposite direction. The corrected RAS-indexed operation exactly
   reproduces Adam's original corrected topology field on all 5,002,140 crop
   voxels when the original WM input is held fixed. It also matches all five
   independently executed selected Java revision-1.7 fixtures exactly. This
   held-input match is component evidence, separate from the independent chain.
   Dura and generic callers retain their existing default paint-mask branch.
2. **Hemisphere separation:** the corrected `run_maCRUISE_v3_2_classern.m` in
   the supplied `MaCRUISE_v3_1_0/` directory enables `SepLeftRight` by default.
   The full MinFilter3D boundary of left/right cortical labels 101/102 is
   excluded from all four masks. All four now exactly match Adam's originals.
3. **TGDM handoff:** the BASE wrapper leaves the optional Brain Mask Level Set
   unset. Python BASE and continuation now do the same. The initial brain mask
   continues to control preprocessing/TOADS. Explicit standalone TGDM bounds
   remain available.
4. **Axial preprocessing:** mask before WINDOWED_SINC, retaining the selected
   table, float coordinate products, double accumulation and disabled clipping.
   The production port exactly matches selected Java on 9,830,400 held source
   voxels. Intensities reach TOADS without another binary remasking step.
5. **Dura order:** run source-sensitive Dura operations in LPS, restoring public
   products to RAS. With original tissue inputs held fixed, mask disagreements
   decrease from 90 voxels to one. The independent full-chain Dura mask improves
   from 166,096 differing voxels to 6,075 (Dice 0.9612702 to 0.9986364).
6. **Historical N4 policy:** explicitly reproduce the old automatic 200-bin
   Otsu estimation mask, and correct the entire image with its estimated bias
   field. The mask matches independent ITK exactly on Adam's scan and thresholds/
   masks match on 40 additional synthetic cases. N4 masked correlation against
   the original rises from 0.9874264 to 0.9990657. The acquisition grid and exact
   stored spacing are retained before working-grid conversion.

Historical N4 policy: [ANTs 2.1.0 source](https://github.com/ANTsX/ANTs/blob/v2.1.0/Examples/N4BiasFieldCorrection.cxx).
The exact original N4/runtime identity is unknown; N4 equality is unproven.
With original second-pass tissue maps held fixed, fusion exactly reproduces
the supplied original GM and WM fields. The independent second-pass TOADS
classification still differs at 248,749 voxels (versus 405,510 in 0.1.1),
locating substantial remaining differences upstream of these matched components.
The precise cause of the remaining independent reconstruction errors is not
established by this evidence.
Standalone N4 and existing PLUS calls retain their current `all` estimation
mask default. This release makes no new claim of PLUS equivalence, validation
on other subjects, or scientific/clinical validity from passing software tests.

Upstream-only 0.1.2 and subsequent partial correction candidates were withheld:
the improved upstream results alone worsened final label agreement. The final
0.1.5 comparison above includes the topology handoff fix and is delivered as an
explicitly investigational testing package, with the remaining regressions
recorded. It must not be presented as a complete fix or overall accuracy upgrade. Patient inputs,
patient-derived images/meshes, original programs, Java/MATLAB/C++ harnesses,
containers, isolated reference dependencies and experimental prototypes are
excluded from this distribution. The installed test log and numerical JSON
receipts are included; local visualization stays in the source checkout.

## Local evidence

- Final chain: `build/adam-upward-native-2026-10-29/run.json`
- Full comparison and visual QA: `build/adam-upward-parity-2026-10-29/`
- Handoff: `build/pymacruise-adam-0.1.5-2026-10-29.tar.gz`
