# TOADS paper-to-code notes

This note connects the published TOADS method to the specific portions of
PyMaCRUISE that currently implement related operations. It is explanatory
documentation, not evidence that PyMaCRUISE implements either paper in full.
For behavior, use the selected VandyCRUISE Java source and its wrapper settings
as the authority. The source selection, wrapper contract, and known deviations
are recorded in [`../toads-port-plan.md`](../toads-port-plan.md) and
[`../toads-atlas-registration.md`](../toads-atlas-registration.md).

## Published papers

1. Pierre-Louis Bazin and Dzung L. Pham, “Topology-preserving tissue
   classification of magnetic resonance brain images,” *IEEE Transactions on
   Medical Imaging*, 26(4), 487–496 (2007).
   [doi:10.1109/TMI.2007.893283](https://doi.org/10.1109/TMI.2007.893283).
   The paper describes TOADS as combining statistical tissue classification,
   topology-preserving fast marching, and image registration to impose
   object-level relationships.
2. Pierre-Louis Bazin and Dzung L. Pham, “Homeomorphic brain image
   segmentation with topological and statistical atlases,” *Medical Image
   Analysis*, 12(5), 616–625 (2008).
   [doi:10.1016/j.media.2008.06.008](https://doi.org/10.1016/j.media.2008.06.008).
   This paper develops the use of topological and statistical atlas priors in
   homeomorphic brain segmentation, with tissue classification and fast
   marching as parts of the framework.

Bibliographic details and abstracts were checked against [PubMed for the 2007
paper](https://pubmed.ncbi.nlm.nih.gov/17427736/) and [the 2008 paper's
publisher/DOI record](https://doi.org/10.1016/j.media.2008.06.008).

## Code crosswalk

| Published-method component | PyMaCRUISE code | What the code contributes | Boundary of the claim |
|---|---|---|---|
| Statistical tissue classification and intensity terms (2007, 2008) | [`toads_preprocessing.py`](../../src/pymacruise/toads_preprocessing.py), [`toads_initialization.py`](../../src/pymacruise/toads_initialization.py), [`toads_math.py`](../../src/pymacruise/toads_math.py) | Source-shaped robust intensity bounds and normalization, initial centroids and exact reciprocal-distance memberships, selected optimized tissue likelihood and membership update kernels. | These are source-anchored primitives and selected passes. They do not by themselves implement or validate the paper's full alternating classification and segmentation procedure. Use the selected Java methods as the equation and edge-case authority. |
| Statistical atlas and topology-template priors (2008) | [`toads_atlas.py`](../../src/pymacruise/toads_atlas.py), [`toads_shapes.py`](../../src/pymacruise/toads_shapes.py), [`toads_crop.py`](../../src/pymacruise/toads_crop.py) | Descriptor/asset loading, rigidly transformed shape-prior sampling, source-convention bounds, crop/pad, and categorical topology-template sampling. | Loading and sampling are not atlas construction, deformable registration, or a whole atlas-driven segmentation. World-origin reconciliation and complete registration remain distinct validation questions. |
| Image registration component (2007, 2008) | [`toads_registration.py`](../../src/pymacruise/toads_registration.py), [`toads_registration_sampling.py`](../../src/pymacruise/toads_registration_sampling.py) | Selected rigid exact-membership objective terms, derivatives/Jacobians, sampled atlas fields, pyramid scheduling, and gradient line-search operations. | Current synthetic source comparisons validate limited formulas and composed update slices. They do not show full-atlas alignment, full-run convergence, or paper-level equivalence. |
| Topology-preserving multiple-object front / fast-marching strategy (2007, 2008) | [`toads_heap.py`](../../src/pymacruise/toads_heap.py), [`toads_thinning.py`](../../src/pymacruise/toads_thinning.py) | The source-compatible priority queue plus selected growing and thinning queue initialization/transitions, including class volume changes and queue scheduling. | The heap and validated events are components, not the full algorithm. A source-matched queue pass on a synthetic phantom does not establish complete VandyCRUISE iteration, convergence, or clinical-data behavior. |
| Local digital-topology and homeomorphic label-change constraint (especially 2008) | [`toads_critical_lut.py`](../../src/pymacruise/toads_critical_lut.py), [`toads_topology.py`](../../src/pymacruise/toads_topology.py) | The selected 18/6 critical-point lookup and local multi-label predicate checking the proposed/old objects and neighboring label relations. | A local gate only preserves the source's local admissibility decision. It is not an independent proof that a complete output volume is homeomorphic to its template. |

`toads_dura.py` implements a separate MaCRUISE inter-pass dura-removal stage
using selected JIST/VandyCRUISE operations. It is not described here as an
implementation of either core TOADS paper. BET and SPECTRE are outside the
selected Python runtime path; the user-provided native-T1 brain mask and
working-space deviation are specified in the TOADS implementation notes.

## Reading the annotations

The source-module docstrings point to this note and name the paper concept
only where code has a defensible relationship to it. “Implements a component”
means that the named operation occupies that role in the published framework;
it does not assert that the source implementation and publication are
identical, that the whole workflow is ported, or that numerical equivalence is
established. Individual validation statements should retain their measured
scope, such as “matches the Java harness reference,” and should cite their
synthetic fixture and full-suite evidence in the corresponding report.

The working-space, atlas, registration, and remaining end-to-end gaps are
tracked in the linked source implementation notes. The selected Java version
and its wrapper are the behavior authority where the papers leave engineering
details unspecified or the implementation has evolved.
