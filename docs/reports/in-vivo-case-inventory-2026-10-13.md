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

The local input folder currently has 41 files (~602 MB at the latest check).
The initial inventory found 22 NIfTIs and read all their headers. The
high-resolution anisotropic FLAIR is 192×512×512 at approximately
1.000005×0.488281×0.488281 mm; its T1 is 176×240×256 at approximately 1 mm
isotropic. A T2, GRE magnitude/phase pair, DWI, BOLD, and field maps are
present. The high-resolution FLAIR needs physical-space registration and
resampling for comparison with T1-grid stages.

The reference directory is spelled `Output/` (capital O), unlike case
100249's `outputs/`. A fresh read-only check on 2026-10-13 found no files in
that directory; only empty `Surfaces_MNI/`, `Surfaces_FreeView/`, and
`SegRefine/` directories were present at the initial inventory. The input
file count and size may change while transfer is ongoing, so recheck before
using this as a complete case. It is currently an additional FLAIR input, not
a saved MaCRUISE reference.

## Validation use

Use 100249's populated files only for T1-chain checkpoint comparisons after
reconciling transforms and stage provenance. Treat both FLAIR images as
registration/robustness inputs until paired source FLAIR intermediates and
transforms are available. Keep all subject data external to the repository;
record only permitted hashes, geometry, comparison metrics, and local paths
in validation reports.
