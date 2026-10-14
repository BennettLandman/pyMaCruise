# MaCRUISE surface-guided segmentation refinement

The selected V3 runner (`run_maCRUISE_v3.m`, SVN r9360) calls
`run_SegmentationRefine.m`, which selects inner/outer CRUISE level sets and
passes the source MaSeg file to `MAsegRefineBySurf.m`. The V3 wrapper sets all
four coefficients to 0.2 mm and enables central processing. The helper
divides each coefficient by the largest MaSeg voxel spacing, uses a 5³
neighborhood (radius 2), and permits up to 50 sequential fill sweeps.

`pymacruise.segmentation_refine.refine_segmentation_by_surfaces` implements
three source rules, in order:

1. Labels `>=100` with `inner < -0.2 mm / max_spacing` are filled with the
   local mode of neighboring labels 44/45.
2. Remaining labels `>=100` with `outer > 0.2 mm / max_spacing` are cleared
   to 0.
3. When central handling is enabled, labels 44/45 beyond the inner threshold,
   plus label-0 voxels between the two thresholds, are filled from the local
   mode of neighboring labels `>=100`.

Fresh execution of the original `MAsegRefineBySurf.m` under MATLAB R2024b Update 5 established that numeric positive-integer `tabulate` includes zero-count rows for every integer from 1 through the largest label. Thus `[~,maxind]=max(foundtable(:,2))` returns the actual modal label for the WM (44/45) and cortical (>=100) candidates. Ties select the lowest label. The earlier claim that the MATLAB source creates labels 1 or 2 was based on treating the table as distinct-values-only and is superseded by this runtime evidence.

The default `source_faithful=True` now assigns the source-equivalent modal label. `source_faithful=False` remains accepted for compatibility and produces the same labels for these candidate groups. Local updates remain sequential in MATLAB column-major index order. The MATLAB helper clips neither its radius-2 boxes nor their bounds, so its edge accesses can fail; Python still clips those boxes. This separate boundary deviation is not validated by the interior reference fixture and must remain explicit.

The API requires co-gridded 3D NIfTI inputs and integer labels. Non-RAS inputs
warn and are reoriented to RAS; the output preserves the segmentation grid
and integer label dtype. The result includes candidate masks, changed masks,
per-rule counts, unresolved candidate counts, sweep counts, and the effective
voxel thresholds so every modification can be audited.

## Published-method context: MaCRUISE CCSE

Huo Y, Plassard AJ, Carass A, Resnick SM, Pham DL, Prince JL, Landman BA.
“Consistent cortical reconstruction and multi-atlas brain segmentation.”
*NeuroImage*. 2016;138:197–210. DOI
[10.1016/j.neuroimage.2016.05.030](https://doi.org/10.1016/j.neuroimage.2016.05.030);
[PubMed](https://pubmed.ncbi.nlm.nih.gov/27184203/);
[open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4927397/).

Section 2.4 describes CCSE as using inner and outer surface level sets plus
consistency distances to reconcile volume labels with reconstructed cortical
surfaces. The selected V3 helper implements that broad relationship: it clears
cortical labels beyond the outer surface, fills the inner-side region from WM
neighbors, and assigns cortical labels between surfaces. Exact thresholds,
label groups, neighborhood size/order, clipping at volume edges, and the
`maxind`/positive-integer `tabulate` behavior are defined by MaCRUISE SVN r9360 and its
V3 wrapper, not inferred from the paper. In particular, the V3 wrapper uses
0.2 mm coefficients here, while the paper's main evaluation uses 0.5 mm.

The MaCRUISE paper also describes multi-atlas label fusion and TOADS-to-label
segmentation fusion (§2.2). The current workflow accepts a MaSeg segmentation
as input; these upstream fusion algorithms are not implemented by this module
and must not be attributed to it. MaACE is a separate cortical-boundary
enhancement method described in §2.3.1 and is tracked under P51.

## Synthetic validation

The committed 11³ `.nii.gz` fixtures in
`tests/data/synthetic/refinement/` use an anisotropic 0.5/0.75/1.0 mm RAS grid.
They cover inner-side WM filling, a tied 44/45 vote, outer-side clearing,
central WM and missing-GM filling, an unresolved volume-border candidate,
and unrelated-label preservation. Tests also cover central processing off,
non-RAS warnings/output, grid mismatch, integer-label validation, the
coefficient conversion on a 0.5 mm grid, and invalid options. No actual MRI
data are used.

Tests cover both compatibility modes and their agreement with the actual MATLAB-written 13³ golden reference in `tests/data/synthetic/matlab_refinement_reference/`. That fixture covers all three refinement rules with interior candidates. Fresh MATLAB probes also cover positive-integer uint16/double labels, missing histogram rows, and tied counts. See [the direct-source comparison](reports/matlab-refinement-direct-source-comparison-2026-10-14.md). Whole-workflow and volume-border comparisons remain open.
