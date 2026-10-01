# ACE and MaACE: paper-to-code guide

This note connects the published method descriptions to the selected PyMaCRUISE
V3 implementation. It distinguishes the algorithmic idea in a paper from
parameters and processing details taken from the supplied VandyCRUISE/JIST
source.

## Papers

1. D. Tosun, M. E. Rettmann, X. Han, X. Tao, C. Xu, S. M. Resnick, D. L. Pham,
   and J. L. Prince, “Cortical surface segmentation and mapping,”
   *NeuroImage*, 23 Suppl. 1 (2004), S108–S118.
   [Publisher DOI: 10.1016/j.neuroimage.2004.07.042](https://doi.org/10.1016/j.neuroimage.2004.07.042)
   · [Author-hosted paper](https://iacl.ece.jhu.edu/~prince/pubs/2004-NI-Tosun-IPAM.pdf).
   ACE is described in the section “Anatomically consistent enhancement,”
   pp. S110–S111, equations (2)–(3).
2. Y. Huo, A. J. Plassard, A. Carass, S. M. Resnick, D. L. Pham, J. L. Prince,
   and B. A. Landman, “Consistent cortical reconstruction and multi-atlas brain
   segmentation,” *NeuroImage*, 138 (2016), 197–210.
   [Publisher DOI: 10.1016/j.neuroimage.2016.05.030](https://doi.org/10.1016/j.neuroimage.2016.05.030)
   · [Open manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC4927397/).
   MaACE is described in §2.3.1, equations (5)–(9).

## Method idea and implementation map

| Paper-level method | PyMaCRUISE location | What the code does and where its exact details come from |
|---|---|---|
| ACE starts a distance field at the GM/WM interface and solves a speed-weighted Eikonal equation. Tosun et al. use `F(x) = 1 - 0.9 μ_CSF(x)` so travel is slower in CSF and the resulting shock skeleton shifts into sulci. | [`ace_update_speed`](../../src/pymacruise/ace.py) and [`_distance_and_characteristic_fmm`](../../src/pymacruise/ace.py) | `ace_update_speed` forms the conventional ACE speed from the CSF field; the FMM calculates the travel-time/distance characteristic. The six-neighbor updates, interface initialization, and float behavior are ported from VandyCRUISE. They are source details, not specifications of the paper. |
| ACE uses the shock/skeleton set to create a thin digital separation between opposing GM banks and suppresses GM membership at those points. | `ace_enhance_arrays` in [`ace.py`](../../src/pymacruise/ace.py) | Gradient magnitude and distance thresholds plus component cleanup create candidate support; the output enhanced-GM image is set to zero on that support. The selected adapter returns a binary CSF/skeleton product; it does not expose a fuzzy CSF map with the removed GM values explicitly accumulated into it as described conceptually in the paper. |
| Huo et al. MaACE generalizes ACE by using the larger of CSF membership and multi-atlas boundary membership in the speed: `F(x) = 1 - 0.9 max(μ_CSF(x), μ_boundary(x))`. The paper defines a cortical boundary from label changes in the 26-neighborhood while excluding neighborhoods containing WM labels. | Separate selected boundary stage: `multi_atlas_csf_boundary` in [`boundary.py`](../../src/pymacruise/boundary.py) and `run_plus_boundary_stage`; precomputed `multi_atlas_csf` / merge placements in `ace_enhance_arrays`, `ace_seed_fields`, and `ace_update_speed` in [`ace.py`](../../src/pymacruise/ace.py); plus ACE wrapper: [`plus_pipeline.py`](../../src/pymacruise/plus_pipeline.py) | **Scope difference:** the selected source has a separate stage that computes a binary `csf_MaBoundary`-style product from atlas labels using V3 `MinFilter3D`; the ACE adapter accepts a precomputed atlas-CSF/boundary image. The Python plus-stage helpers do not themselves establish an end-to-end handoff between those functions. The relation of that product to Huo et al.'s exact Eq. (7) `μ_boundary` is not established. The base adapter promotes positive input voxels to CSF=1 after diffusion; the plus adapter promotes positive voxels to speed=0.1 after the speed Gaussian. These selected-wrapper merges are separately compared with Java, but they do not apply or validate Eq. (7)'s max-based speed rule using the paper-defined boundary map. Do not describe this path as a complete MaACE implementation. |
| The separating support is thinned to one voxel. | [`thin_ace`](../../src/pymacruise/skeleton.py) | The need for a thin skeleton is in the ACE/MaACE descriptions. The actual distance-priority peeling, simple-point and endpoint predicates, connectivity semantics, iteration order, and 254 encoding are traced to the supplied VandyCRUISE `ThinACE` source. This cleanup is distinct from MaACE's boundary-prior calculation. |

## Source-specific choices in this selected V3 path

The following are not asserted to be prescribed by either paper. They are
parameters or operation order taken from the selected V3 adapter and
VandyCRUISE algorithm classes:

- The input interface is initialized from `intensity_threshold - white`, with
  a default threshold of 126.9. The papers describe a 0.5 isosurface of WM
  membership; this API receives source WM intensity values, not a normalized
  WM membership image.
- Initial CSF seeding uses WM `>= 30`, combined WM+GM `>= 210`, a four-iteration
  18-neighbor closing, and `sqrt(naive_csf / 254)` on eligible voxels.
- The adapter diffuses initialized CSF, merges positive atlas-CSF values after
  diffusion, and applies a recursive Gaussian. The chosen alphas and iteration
  counts are recorded in [`ace-port.md`](../ace-port.md).
- Narrow-band evolution, weighted-distance cutoffs, gradient cutoff, component
  size cutoff, and border clearing follow the traced V3 source configuration.
- Normal tensor diffusion preserves a finite initialized eigenvector when a
  zero tensor is encountered. `source_faithful=True` is available on the
  stage-level diffusion function to reproduce the selected Java code's NaNs.
  This is an explicitly documented stability deviation, not part of either
  paper's ACE/MaACE definition.

The selected Java wrapper's stage-by-stage comparison is documented in
[`ace-port.md`](../ace-port.md). That evidence establishes agreement with that
synthetic Java harness reference for the tested path; it does not establish
equivalence to every publication equation or input configuration.
