# In-vivo case inventory — 2026-10-13

This inventory records available local data for later real-image validation.
It is not a numerical comparison, anatomical acceptance result, or proof that
each folder is a complete run. No images or subject outputs are copied into
PyMaCRUISE or committed.

## Case 100249 (`sub-010249`)

The local case contains 109 files (~1.11 GB): 37 source/auxiliary files under
`/Users/landmaba/dev/mipav/sample-macruise/100249/sub-010249/ses-01` and 72
saved files under `outputs/`. The initial inventory read 69 NIfTI headers and
found no zero-length files.

Its FLAIR is a low-resolution anisotropic image (384×512×28;
0.4296875×0.4296875×5 mm). The MP2RAGE T1 is 176×240×256 at 1 mm isotropic;
T2 is 176×256×256 at 1 mm isotropic. DWI, BOLD, and field-map data are also
present.

The populated saved output directory contains candidate T1-chain checkpoints:
N4, MultiAtlas stripping, two TOADS passes, Dura removal, cortical separation,
TGDM inputs, refinement, surfaces, thickness, and MNI/freeview derivatives.
The output name and saved artifacts establish a T1 MaCRUISE run; they do not
establish a FLAIR/PLUS reference. No FLAIR-to-T1 transform, registered
bias-corrected FLAIR, lesion mask, or run configuration log was identified.
Volume dimensions alone are insufficient for voxel comparison: qform/sform
and physical transforms vary across saved stages.

## Case 100260 (`sub-010260`)

The local case now has 114 files (~1.20 GB), including 41 BIDS input files
and 71 NIfTIs. The high-resolution anisotropic FLAIR is 192×512×512 at
approximately 1.000005×0.488281×0.488281 mm; its T1 is 176×240×256 at
approximately 1 mm isotropic. A T2, GRE magnitude/phase pair, DWI, BOLD, and
field maps are present. The high-resolution FLAIR needs physical-space
registration and resampling for comparison with T1-grid stages.

A same-grid BET mask and a SLANT-TICV segmentation are available. Independent
header/value checks confirm that the T1, mask, and SLANT labels are finite,
RAS-oriented, and share the same 176×240×256 grid; the mask is binary and the
SLANT map contains integer BrainCOLOR labels. These are usable inputs for a
Python BASE run.

**Reference update (2026-10-13):** the previously observed `Output/` folder
containing a copy of case 100249's outputs is no longer present. The saved
`MaCRUISEv3.2.0/logs/macruise.log` records MultiAtlas skull stripping ending,
then a MaCRUISE wrapper exit status 1; required `SEG_refine.nii.gz` and the
thickness MAT input to the PDF stage were missing. No completed, subject-
specific source output was found for 100260, so do not score it as an
independent source reference. A fresh Python BASE run using this subject's T1,
BET mask, and SLANT labels is underway in ignored `build/`; it is execution
evidence only, not source equivalence. Its output and anatomy audit must be
reviewed after completion. The 100260 FLAIR remains an input for later PLUS
validation, not an existing PLUS reference.

## Validation use

Use 100249's populated files only for T1-chain checkpoint comparisons after
reconciling transforms and stage provenance. Treat both FLAIR images as
registration/robustness inputs until paired source FLAIR intermediates and
transforms are available. Keep all subject data external to the repository;
record only permitted hashes, geometry, comparison metrics, and local paths
in validation reports.
