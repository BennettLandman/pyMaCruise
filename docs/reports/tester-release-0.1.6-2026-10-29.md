# PyMaCRUISE BASE tester release 0.1.6 — October 29, 2026

Focused diagnostic/engineering testing only; an overall accuracy improvement over 0.1.1 is not established.

The candidate completed a fresh full raw-input run in a clean environment on Adam's supplied case, and all output/grid checks passed. Results remain mixed: refined-label disagreements are 115,560 versus 97,837 in 0.1.1, and central mean surface distance is 0.4069 versus 0.3465 mm. This is a focused diagnostic/engineering testing release, not a completed equivalence result or an overall accuracy upgrade. Please retain the previous outputs for comparison.

## Engineering execution evidence

The delivered wheel was installed in a separate clean Python 3.13.17 environment
with its declared dependencies. All **128 targeted installed-wheel checks pass**
(17 warnings), and **eight tester-tool checks pass**. These selections are
separate; the entire repository suite was not rerun. Dependency checking reports
no broken requirements. The tested scientific versions are pinned in
`constraints-tested.txt` and recorded in `release.json`.

This release completed a **fresh, monolithic raw-input BASE run** on the supplied
T1, brain mask and MaSeg, using this wheel and its full default settings. No
original scientific intermediate was supplied to that reconstruction and no
upstream checkpoint from an earlier Python version was substituted. All ten
stages completed. TOADS converged naturally at iterations
19 and 15.
The packaged launcher then successfully resumed that completed run, verified
and skipped its stages, checked all 84 stage artifacts, and created feedback.
The same installed package produced the final original-output comparison.

All five final NIfTIs preserve the original T1 shape, affine, qform/sform and
storage axes. SEG_refine uses integer storage with no fractional labels.
Primary meshes use physical RAS millimetres; explicitly named voxel-index
comparison copies preserve anatomy through the known reference affine. The
feedback tool was checked on both failed/incomplete and completed runs; its
archive excludes volumes, meshes, arrays, screenshots and external symlinks.
It records package/platform details and the ANTs native-library identity.

The written voxel-coordinate copies retain every mesh face and round-trip to
the primary physical coordinates within 5.02e-10 mm. This macOS runtime emits
NumPy floating-point warnings in the affine matrix multiplication used by the
comparison tools; an independent explicit non-BLAS evaluation confirms finite
coordinates agreeing within 2.85e-14 mm. The numeric receipts are included as
`validation/voxel-copy-qa.json` and `validation/affine-arithmetic-qa.json`.

These are engineering/completion checks, separate from numerical agreement.

## Same-case numerical comparison

The following results use the same supplied inputs and original MaCRUISE
reference. Original raw VTKs are mapped through Raw/Seg's affine and verified
against the original ASC export; no fitted registration, smoothing, component
filtering or cosmetic label cleanup is used. 0.1.1 is the earlier full-default,
native-grid release. 0.1.5 includes the forced historical N4 policy and other
source fixes; a separate clean-environment 0.1.5 control reproduced its previous
results. 0.1.6 restores the current unmasked N4 default while retaining the
other source fixes.

| Label measure | 0.1.1 | 0.1.5 | 0.1.6 |
|---|---:|---:|---:|
| Differing refined-label voxels | 97,837 | 125,524 | 115,560 |
| Labels below Dice 0.97 | 84 | 61 | 72 |
| Mean nonzero-label Dice | 0.93822 | 0.93193 | 0.92689 |

| Surface | Mean mm 0.1.1 | Mean mm 0.1.5 | Mean mm 0.1.6 | 0.1.6 sampled max mm | 0.1.6 components |
|---|---:|---:|---:|---:|---:|
| Inner | 0.36073 | 0.55236 | 0.52078 | 10.19276 | 1 |
| Central | 0.34652 | 0.42910 | 0.40685 | 9.74016 | 1 |
| Outer | 0.26585 | 0.17433 | 0.18718 | 9.84090 | 1 |

Means average both directed all-vertex-to-opposite-triangle means. The sampled
maximum is not a continuous triangle-to-triangle Hausdorff bound. Original
meshes have one component each. Project numerical gates require every nonzero
label Dice ≥0.97, symmetric mean surface distance ≤0.01 mm and sampled maximum
≤1 mm. The combined result is **False**. Detailed per-label measurements
are retained in `validation/comparison.json`; passing execution checks does
not establish source, scientific or clinical equivalence.

## Changes and controlled component evidence

The full default schedule remains registration 4/50, post-refresh/max
registration 5/5, TOADS up to 50, Dura 5, ACE 25, GVF 200 and TGDM 6/4/7,
with WM iso-level 0.5, ACE boundary mode 2 and source-faithful refinement.
Native export restores the original acquisition grid. There are no short
iteration overrides in the supplied launcher.

N4 no longer silently forces a historical automatic Otsu mask. It defaults to
`all`, following the current installed unmasked backend. The source wrapper
supplies no explicit mask, and its exact original runtime remains unknown.
`--n4-estimation-mask legacy-otsu` permits an explicit historical-policy
comparison in a separate fresh directory. The policy is recorded and checked
in run signatures. The default N4 output exactly matches 0.1.1's current bare
backend output at all 9,830,400 voxels. A closer intermediate N4 correlation
alone was insufficient evidence to prefer the previous forced policy.

The source-derived axial interpolation, LPS Dura processing, v3.2 hemisphere
boundary exclusion, upward topology correction and null optional TGDM bound
remain enabled. Earlier held-input tests exactly reproduce all four cortical
masks, the GM/WM fusion products and the 5,002,140-voxel corrected topology crop.
Additional frozen-0.1.5 held-input diagnostics now exactly reproduce the
original ACE enhanced GM and thinned skeleton on that crop. These downstream
operator implementations are unchanged in 0.1.6; their held-input evidence is
separate from its independent raw-input reconstruction.

The original reference-space GVF export retained only its first vector
component and omitted inverse-crop high planes/out-of-image padding. On the
voxels actually retained by that export, Python's first component has mean
absolute error 6.93e-8 and maximum 8.64e-7. Full original three-component equality
cannot be established from that scalar export. With the original tissue/ACE/
skeleton inputs and reconstructed Python GVF, inner and central zero-level
signs match exactly; outer sign Dice is 0.9999904 with 29 disagreements. These
are controlled diagnostics, not a fresh original full-workflow replay. Raw
scalar differences outside the source-export coverage are not anatomical
surface error measurements.

The clean-environment end-to-end case, output checks, original-output comparison,
launcher/resume and feedback tests support the stated tester scope. Broader
subjects/platforms and exact original-runtime equivalence remain unverified.

## Tester handoff

Use the bundled README, `run-base.py`, `compare-with-macruise.py` and
`collect-feedback.py`. The email draft contains an explicit download-link
placeholder for the project owner. Nothing was emailed or uploaded. The
package contains numeric validation receipts and notices, with no patient
images/meshes or original Java/MATLAB/C++ programs, containers or reference
harnesses. External Python dependencies are installed through the wheel's
`ants` extra under the supplied tested constraints.

## Local artifacts

- Full fresh run: `build/adam-tester-final-0.1.6-2026-10-29/run.json`
- Release evidence: `build/adam-release-validation-2026-10-29/`
- Archive: `build/pymacruise-testers-0.1.6-2026-10-29.tar.gz`
- Email: `docs/testing/tester-email-0.1.6.md`
