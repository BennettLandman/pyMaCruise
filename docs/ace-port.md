# ACE port progress

This document tracks the Python implementation of the selected MaCRUISE V3
Anatomically Consistent Enhancement (ACE) path. The parameter contract is read
from `run_Ace_jist.m`, `MaACE.java`, and
`MaAnatomicallyConsistentEnhance.java`; algorithm details are read from the
VandyCRUISE ACE classes. Their notices and attribution inventory remain under
`licenses/upstream/` and in `license_info.md`.

## Selected source configuration

- WM threshold for the initial level set: 126.9.
- Gaussian alpha: 1.8.
- Structure-tensor/diffusion alpha: 1.6, 25 diffusion iterations.
- Updated multi-atlas CSF: merge positive values after diffusion.
- ACE narrow band: time step 0.1, one outer iteration, 10 inner updates,
  curvature threshold 5.0, narrow band `|phi| <= 6`, reinitialization range 4.
- Final membership: exclude weighted distance `<= 0.25` and prohibited voxels;
  use gradient magnitude threshold 0.85, then remove 18-connected components
  of size `<= 60`, promote remaining memberships to 254, and clear the image
  border from both outputs.

## Implemented

`src/pymacruise/ace.py` currently provides:

1. ACE input validation and binary WM closing with the source 18-neighbor
   footprint and border behavior.
2. Initial CSF speed and prohibited-region construction using source thresholds.
3. The Deriche recursive Gaussian and the structure-tensor anisotropic
   diffusion update, including source boundary differences and ten power-method
   iterations. Zero tensors retain the initialized basis vector to avoid the
   source implementation's divide-by-zero NaNs.
4. The after-diffusion positive multi-atlas-CSF merge and speed-function blur.
5. Unit-speed signed fast marching with subvoxel interface initialization,
   source six-neighbor Eikonal updates, and distance caps.
6. Narrow-band ACE curvature evolution with synchronous updates and signed
   distance reinitialization.
7. Speed-weighted fast marching for travel-time and distance-characteristic
   fields, final ACE membership thresholds, 18-connected cleanup, enhanced-GM
   suppression, and border clearing.
8. Array and NIfTI APIs. The NIfTI API validates a common physical grid,
   reorients non-RAS inputs with a warning, and returns RAS NIfTI images.
   The skeleton returned here is the cleaned pre-thinning ACE skeleton;
   thinning is a separate roadmap step (P27).

## Synthetic evidence

`tests/test_ace.py` has 13 deterministic synthetic checks covering the Gaussian
filter, morphology, component cleanup, thresholds, diffusion invariants,
source boundary differences, multi-atlas merge placement, subvoxel planar
distance, narrow-band plane evolution, weighted planar travel time, complete
array-pipeline repeatability, and RAS NIfTI output/warnings. The complete
offline suite currently passes 83 tests. No real imaging scans have been used
for this port.

## Remaining validation before P26 can pass

- Generate an independent trusted reference from the source Java implementation
  (or a hand-audited minimal Java harness) for every numerical ACE stage; compare
  arrays and report MAE, maximum absolute error, Dice, and surface distances.
- Confirm `DistanceFieldAndSkeleton` tie handling and propagated
  distance-characteristic values against that reference, especially where the
  speed is zero or several interface axes cross at one voxel.
- Compare all required segmentation regions against the project acceptance
  thresholds: each Dice score at least 0.97, symmetric Hausdorff distance at
  most 1 mm, and symmetric mean nearest-surface distance at most 0.01 mm.
- Exercise the default 25 diffusion iterations and default atlas merge on
  larger synthetic volumes; characterize memory and runtime before any actual
  scan is introduced.

Consequently, implementation progress is substantial but P26 remains **Started**
until source-reference numerical evidence and all acceptance criteria pass.
