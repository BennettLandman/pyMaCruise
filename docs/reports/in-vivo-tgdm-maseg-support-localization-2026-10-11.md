# In-vivo TGDM support localization — 2026-10-11

## Result

The human source-profile TGDM terminal surfaces cover **81.963%** of the
selected input MaSeg GM and **83.253%** of WM. This is a support measure, not
a Dice score. The prescribed quality threshold is at least 97% for each
tissue.

The TGDM crop contains all selected MaSeg GM and WM voxels, so crop truncation
accounts for none of the gap. The T1-grid user mask excludes only 23,511 of
642,291 GM voxels (3.66%) and 132 of 428,303 WM voxels (0.031%). Most
unsupported tissue is inside the supplied mask and outside the relevant
surface region:

| MaSeg support partition | Voxels | Interpretation |
|---|---:|---|
| GM inside inner surface | 11,203 | Beyond the inner boundary |
| GM between inner and outer surfaces | 526,441 | Supported cortical ribbon |
| GM outside outer surface | 104,647 | Outside the outer boundary |
| WM inside inner surface | 356,577 | Supported white matter |
| In-mask WM outside inner surface | 71,594 | White matter not enclosed by the inner boundary |
| WM outside supplied mask | 132 | Excluded by user-provided mask |

The terminal surfaces are nested at every selected GM voxel; no crossed-sign
GM voxels explain the gap. Of 115,850 unsupported GM voxels, 92,398 are inside
the supplied mask. Of 71,726 unsupported WM voxels, 71,594 are inside the
mask. The final ribbon also includes 90,041 in-mask voxels whose labels are
outside the selected MaSeg GM set. These counts point to a shell-placement or
upstream tissue-label/transfer mismatch, rather than crop truncation or mask
exclusion.

The deficits are distributed across both hemispheres. Unsupported GM spans
approximately RAS bounds `[-65.5,-78.5,-11.5]` to `[64.5,81.5,107.5]` mm;
unsupported WM spans `[-60.5,-76.5,-6.5]` to `[59.5,76.5,102.5]` mm. The
largest WM deficits are label 45 (36,294 voxels) and label 44 (33,016); the
next largest are labels 61 (950) and 62 (928). The largest GM deficits are
labels 142 (4,251), 143 (4,130), 177 (3,714), 190 (3,248), 182 (3,183), 198
(3,099), and 176 (3,011).

## Interpretation and next comparison

The selected MaCRUISE-specific Java `MaGenericTGDM` output was measured against
the same fixed MaSeg labels. Java supports **526,432/642,291 GM (81.9616%)**
and **356,577/428,303 WM (83.2534%)**. Python supports 526,441 GM (81.9630%)
and the same 356,577 WM. Inner sign masks match exactly; outer and derived
ribbon sign masks differ at only 41 voxels (Dice 0.99997923 and 0.99996675).
The Java TGDM input arrays in this comparison were the saved Python-prepared
inputs, so the close surface match isolates TGDM behavior but does not certify
the upstream MaSeg/TOADS preparation.

Source review found that V3 derives GM, WM, mask, and brain-tissue maps from
MaSeg, nearest-neighbor transfers second-pass TOADS maps to the MaSeg/T1 grid,
uses MaSeg GM for the shared 5-voxel-border crop, then calls the existing
`CombineMAandTOADS` stage before topology cleanup and surfaces. The current
Python base chain skips that fusion and crops directly from TOADS GM. The
MaSeg label sets themselves match the MATLAB source. Thus the remaining low
support is upstream of TGDM and is plausibly explained by this missing
preparation/fusion stage. The next step is to compose those already-ported
operations, validate the fused maps against Java, and then rerun only the
TGDM stage on the corrected human inputs.

This is a diagnostic of one supplied human scan, not a change to source
parameters and not a clinical-validity result. All human data and detailed
arrays remain under the ignored `build/invivo-base-2026-10-10-ras-grid/`
directory.
