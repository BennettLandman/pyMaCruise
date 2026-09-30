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
| Statistical tissue classification and intensity terms (2007, 2008) | [`toads_preprocessing.py`](../../src/pymacruise/toads_preprocessing.py), [`toads_initialization.py`](../../src/pymacruise/toads_initialization.py), [`toads_lesion_centroid_initialization.py`](../../src/pymacruise/toads_lesion_centroid_initialization.py), [`toads_math.py`](../../src/pymacruise/toads_math.py), [`toads_exact_dispatch.py`](../../src/pymacruise/toads_exact_dispatch.py), [`toads_default_atlas_membership.py`](../../src/pymacruise/toads_default_atlas_membership.py), [`toads_default_approximate_initial_membership.py`](../../src/pymacruise/toads_default_approximate_initial_membership.py), [`toads_default_approximate_membership.py`](../../src/pymacruise/toads_default_approximate_membership.py), [`toads_atlas_type_membership.py`](../../src/pymacruise/toads_atlas_type_membership.py), [`toads_default_outlier_membership.py`](../../src/pymacruise/toads_default_outlier_membership.py), [`toads_optimized_dispatch.py`](../../src/pymacruise/toads_optimized_dispatch.py), [`toads_optimized_outmask_membership.py`](../../src/pymacruise/toads_optimized_outmask_membership.py), [`toads_lesion_membership_pass.py`](../../src/pymacruise/toads_lesion_membership_pass.py), [`toads_outmask_membership.py`](../../src/pymacruise/toads_outmask_membership.py), [`toads_object_centroids.py`](../../src/pymacruise/toads_object_centroids.py), [`toads_object_variance.py`](../../src/pymacruise/toads_object_variance.py), [`toads_object_image_variance.py`](../../src/pymacruise/toads_object_image_variance.py), [`toads_pair_similarity.py`](../../src/pymacruise/toads_pair_similarity.py), [`toads_rayleigh.py`](../../src/pymacruise/toads_rayleigh.py), [`toads_rayleigh_energy.py`](../../src/pymacruise/toads_rayleigh_energy.py), [`toads_rayleigh_membership.py`](../../src/pymacruise/toads_rayleigh_membership.py), [`toads_inhomogeneity.py`](../../src/pymacruise/toads_inhomogeneity.py), [`toads_inhomogeneity_fit.py`](../../src/pymacruise/toads_inhomogeneity_fit.py), [`toads_inhomogeneity_state.py`](../../src/pymacruise/toads_inhomogeneity_state.py), [`toads_initial_fielded_membership.py`](../../src/pymacruise/toads_initial_fielded_membership.py) | Source-shaped robust intensity bounds and normalization, initial centroids and reciprocal-distance memberships, optional top-k seeds and atlas-type dispatch, selected OBJECT/LESION/OUTMASK/OUTLIER/Rayleigh terms, a bias-field basis and bounded degree-one field fit, Java-style per-channel field attachment state, and a fixed-state composition that attaches fields before the selected initial 11-class membership pass. | The publications motivate statistical tissue classification and atlas-guided segmentation; they do not specify this selected solver's exact field polynomial, setter semantics, or startup call sequence. The degree-one fit, reference-preserving per-channel assignment plus `useField` flag, and membership-after-attachment composition follow selected VandyCRUISE source. Lesion/OUTMASK passes and the Rayleigh distance kernel are similarly source-specific. These bounded operations do not estimate every statistic or implement the papers' complete alternating procedure. Use the selected Java methods as equation and edge-case authority. |
| Selected fixed-degree Chebyshev bias-field fit (P18cya, in progress) | [`toads_inhomogeneity.py`](../../src/pymacruise/toads_inhomogeneity.py), [`toads_inhomogeneity_fit.py`](../../src/pymacruise/toads_inhomogeneity_fit.py) | The selected `InhomogeneityCorrection2.computeApprox3DPolynomialField(n)` source evaluates a 3-D Chebyshev field basis and fits the effective OPTIMIZED_IMAGE/GLOBAL/CHEBYSHEV/FIXED/EQUAL_VARIANCE profile with q=2. Degree one is the current validated fit; P18cya extends the same fixed-degree fit to degrees two and three, with 10 and 20 source-ordered coefficients. | The 2007/2008 papers provide broad statistical tissue-classification/intensity-model context; this precise basis, constructor switch, every-third-voxel fitting scan, p=1 PowerTable weights, class-0 candidate exclusion, system normalization, solve, full-grid evaluation, mask mean-centering, and unit field outside the mask are selected VandyCRUISE implementation details. P18cya's validation gate is a pinned Java method harness on a folded anisotropic RAS `.nii.gz` phantom comparing degree-specific normal systems, coefficients, evaluated fields, centering, and outside-mask units. Until that evidence is recorded, no degree-two/three parity claim is made. The fixed-degree fits do not port the polynomial-degree ramp, field application, multiple-channel correction orchestration, or full TOADS correction workflow. |
| Source-specific startup composition (wrapper lifecycle; not a separate paper claim) | [`toads_initial_prepared_membership.py`](../../src/pymacruise/toads_initial_prepared_membership.py), with [`toads_initial_relation_state.py`](../../src/pymacruise/toads_initial_relation_state.py) | Wires one field attachment per channel → relation, boundary-distance, and closest-label preparation from constructor-zero boundary state → the fixed-input 11-class exact-membership pass. | Variances, image variances, centroid similarities, transformed class-shape maps, and atlas priors are supplied as fixed upstream values; this adapter does not estimate/sample them or iterate TOADS. Its composition test checks ordering and output wiring over already Java-reference-validated components; it is not a new combined Java oracle or evidence of full-workflow equivalence. |
| Statistical atlas and topology-template priors (2008) | [`toads_atlas.py`](../../src/pymacruise/toads_atlas.py), [`toads_shapes.py`](../../src/pymacruise/toads_shapes.py), [`toads_crop.py`](../../src/pymacruise/toads_crop.py) | Descriptor/asset loading, rigidly transformed shape-prior sampling, source-convention bounds, crop/pad, and categorical topology-template sampling. | Loading and sampling are not atlas construction, deformable registration, or a whole atlas-driven segmentation. World-origin reconciliation and complete registration remain distinct validation questions. |
| Image registration component (2007, 2008) | [`toads_registration.py`](../../src/pymacruise/toads_registration.py), [`toads_registration_sampling.py`](../../src/pymacruise/toads_registration_sampling.py), [`toads_approximate_registration.py`](../../src/pymacruise/toads_approximate_registration.py), [`toads_default_approximate_registration.py`](../../src/pymacruise/toads_default_approximate_registration.py), [`toads_default_single_scale_registration.py`](../../src/pymacruise/toads_default_single_scale_registration.py), [`toads_default_registration_handoff.py`](../../src/pymacruise/toads_default_registration_handoff.py), [`toads_lesion_atlas_registration_handoff.py`](../../src/pymacruise/toads_lesion_atlas_registration_handoff.py), [`toads_post_registration_template_refresh.py`](../../src/pymacruise/toads_post_registration_template_refresh.py) | Selected rigid exact-membership objective terms, derivatives/Jacobians, sampled atlas fields, initial approximate and later scale-one registration handoffs, pyramid scheduling, line-search operations, membership-to-registration handoffs, and the source wrapper's post-registration recrop/template refresh. The 13-class handoff includes a LesionTOADS atlas descriptor through the selected base-TOADS registration dispatch. | These bounded compositions connect to the papers' atlas-registration role, but pyramid and wrapper lifecycle details are selected VandyCRUISE behavior. They are not complete atlas alignment or evidence that all classes/wrappers behave equivalently. Synthetic comparisons validate only the reported class order, selected update slices, and handoff state. |
| Topology-preserving multiple-object front / fast-marching strategy (2007, 2008) | [`toads_heap.py`](../../src/pymacruise/toads_heap.py), [`toads_thinning.py`](../../src/pymacruise/toads_thinning.py), [`toads_relations.py`](../../src/pymacruise/toads_relations.py), [`toads_control.py`](../../src/pymacruise/toads_control.py), [`toads_post_registration.py`](../../src/pymacruise/toads_post_registration.py), [`toads_next_default_membership.py`](../../src/pymacruise/toads_next_default_membership.py), [`toads_iteration.py`](../../src/pymacruise/toads_iteration.py), [`toads_iteration_slice.py`](../../src/pymacruise/toads_iteration_slice.py) | The source-compatible priority queue, selected relation/closest-label refresh, growing/thinning events, post-registration stop/continue branch, and a fixed-statistics continuing pass into the next default membership sweep. | These components occupy the paper's broad iterative/topology-aware role, but exact loop thresholds and stage order follow the selected wrapper/source. They are bounded slices: statistics estimation, all wrapper options, repeated registration/iteration to convergence, and the complete published workflow are not established. Java-reference claims remain limited to each recorded synthetic fixture and method boundary. |
| Local digital-topology and homeomorphic label-change constraint (especially 2008) | [`toads_critical_lut.py`](../../src/pymacruise/toads_critical_lut.py), [`toads_topology.py`](../../src/pymacruise/toads_topology.py) | The selected 18/6 critical-point lookup and local multi-label predicate checking the proposed/old objects and neighboring label relations. | A local gate only preserves the source's local admissibility decision. It is not an independent proof that a complete output volume is homeomorphic to its template. |

The selected lesion-bearing atlas path is split across
[`toads_lesion_atlas_type_membership.py`](../../src/pymacruise/toads_lesion_atlas_type_membership.py)
and [`toads_lesion_atlas_registration_handoff.py`](../../src/pymacruise/toads_lesion_atlas_registration_handoff.py).
The first applies the selected base-TOADS 13-class ordering and its LESION
membership branch; the second feeds those fixed-input memberships into the
base-TOADS rigid registration path. This composition does not make the
LesionTOADS iterative engine part of base TOADS or implement full TOADS.

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
