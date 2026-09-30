# Lesion-TOADS paper-to-code notes

## Published method

Shiee N, Bazin P-L, Ozturk A, Reich DS, Calabresi PA, Pham DL. “A topology-
preserving approach to the segmentation of brain images with multiple
sclerosis lesions.” *NeuroImage* 49(2), 1524–1535 (2010).
[doi:10.1016/j.neuroimage.2009.09.005](https://doi.org/10.1016/j.neuroimage.2009.09.005)
· [PubMed](https://pubmed.ncbi.nlm.nih.gov/19766196/) ·
[open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2806481/).

The paper extends TOADS to multichannel images with an additional lesion
class, while preserving the topology of the combined white-matter/lesion
region. It uses topological and statistical atlases and describes modified
relationship terms and channel weighting for lesion-aware segmentation. Those
are descriptions of the published algorithm, not claims that the current
Python modules implement the complete algorithm.

## Implemented-code crosswalk

| Python file | Implemented portion | Connection to the paper | Boundary / validation |
| --- | --- | --- | --- |
| Atlas, constructor state, and preprocessing: [`lesion_atlas_priors.py`](../../src/pymacruise/lesion_atlas_priors.py), [`lesion_state.py`](../../src/pymacruise/lesion_state.py), [`lesion_toads_preprocessing.py`](../../src/pymacruise/lesion_toads_preprocessing.py), [`lesion_covariance.py`](../../src/pymacruise/lesion_covariance.py) | Parse descriptor priors, initialize selected Java constructor buffers/statistics, normalize multichannel inputs, initialize lesion/black-hole centroids, and expose the selected covariance update. | These are statistical-atlas and initialization pieces for the multichannel lesion model. | The descriptor parser does not load topology/shape volumes. Covariance corresponds to a dormant Java method, not an active wrapper stage. P39a/P39b and the covariance slice compare only their named source operations; zero robust spans intentionally raise instead of reproducing Java NaN/infinity. |
| Similarity, statistics, and data terms: [`lesion_similarity.py`](../../src/pymacruise/lesion_similarity.py), [`lesion_statistics.py`](../../src/pymacruise/lesion_statistics.py), [`lesion_statistics_refresh.py`](../../src/pymacruise/lesion_statistics_refresh.py), [`lesion_data_terms.py`](../../src/pymacruise/lesion_data_terms.py), [`lesion_data_energy_maps.py`](../../src/pymacruise/lesion_data_energy_maps.py) | Calculate selected class/lesion similarity, centroid/variance/channel-weight refreshes, and ordinary/lesion/black-hole data-energy terms. P39v adds the selected field-multiplied image path for those data terms. | Supplies intensity and parameter updates for the paper's statistical classification model. | These are separate fixed-input or method-level slices. Special FLAIR/Sulcal-CSF branches, estimation dependencies, and any unsupplied wrapper option must be called out at the corresponding row; they do not establish the whole paper energy model. |
| Membership and objective updates: [`lesion_membership.py`](../../src/pymacruise/lesion_membership.py), [`lesion_gm_ventricle_membership.py`](../../src/pymacruise/lesion_gm_ventricle_membership.py), [`lesion_membership_pass.py`](../../src/pymacruise/lesion_membership_pass.py), [`lesion_objective.py`](../../src/pymacruise/lesion_objective.py), [`lesion_objective_membership_pass.py`](../../src/pymacruise/lesion_objective_membership_pass.py), [`lesion_relation_aware_membership.py`](../../src/pymacruise/lesion_relation_aware_membership.py), [`lesion_statistics_membership_bridge.py`](../../src/pymacruise/lesion_statistics_membership_bridge.py) | Compose selected smoothing, priors, relation gating, distance modulation, OUTMASK handling, normalization, objective reduction, and fixed-input statistic-to-membership updates. | These terms implement parts of the paper's lesion-aware class competition and auxiliary lesion/black-hole map updates. | Each adapter accepts fixed maps or statistics at its documented boundary; unsupported class branches, skipped-shape behavior, and objective exclusions remain local limits. A composition of validated kernels is not the complete iterative algorithm. |
| Atlas registration: [`lesion_shapes.py`](../../src/pymacruise/lesion_shapes.py), [`lesion_registration_handoff.py`](../../src/pymacruise/lesion_registration_handoff.py) | Sample lesion-bearing class shape priors and compose selected updated memberships into atlas registration. | Supplies topological/statistical-atlas alignment components described by the paper. | These source-shaped registration slices do not establish full registration over all deployed atlas assets, wrapper variants, or complete segmentation runs. |
| Relations, topology propagation, and iteration: [`lesion_relations.py`](../../src/pymacruise/lesion_relations.py), [`lesion_growth.py`](../../src/pymacruise/lesion_growth.py), [`lesion_distance.py`](../../src/pymacruise/lesion_distance.py), [`lesion_brainstem.py`](../../src/pymacruise/lesion_brainstem.py), [`lesion_csf_distance.py`](../../src/pymacruise/lesion_csf_distance.py), [`lesion_gm_distance.py`](../../src/pymacruise/lesion_gm_distance.py), [`lesion_ventricle.py`](../../src/pymacruise/lesion_ventricle.py), [`lesion_post_registration.py`](../../src/pymacruise/lesion_post_registration.py), [`lesion_continuing_iteration.py`](../../src/pymacruise/lesion_continuing_iteration.py) | Build selected relation lists and diffusion, run fixed-input growing/thinning operations and capped distance/name-selection helpers, then compose limited post-registration and continuing-iteration slices. | These components correspond to the paper's topology-preserving competition and lesion-class relationship roles. | The capped six-connected distance is a source helper, not the paper's fast-marching evolution. The compositions cover only their documented branches and fixed inputs; they do not show full-volume topology preservation or workflow convergence. |
| Classification and wrapper outputs: [`lesion_classification_update.py`](../../src/pymacruise/lesion_classification_update.py), [`lesion_membership_classification.py`](../../src/pymacruise/lesion_membership_classification.py), [`lesion_subtype_map.py`](../../src/pymacruise/lesion_subtype_map.py), [`lesion_cruise_inputs.py`](../../src/pymacruise/lesion_cruise_inputs.py), [`lesion_dura_outputs.py`](../../src/pymacruise/lesion_dura_outputs.py) | Transfer selected class/lesion/subtype labels and prepare the selected CRUISE/dura-removal output arrays. | Classification supports lesion-aware reporting in the algorithm family; output channel packing is a wrapper/data-interface responsibility. | Output packing and dormant classification helpers are not paper-method evidence. The selected wrappers do not call every ported helper; active reachability is documented separately. |

## Version and completeness limits

The project selected VandyCRUISE as the source tree for the port, but the
exact deployed plugin version used by historical MaCRUISEplus runs is still
unknown. Consequently these mappings describe the cited 2010 method at a
conceptual level and the specific methods/data parsed from the selected
VandyCRUISE source; they do not establish that the selected tree is identical
to the paper implementation or the historical deployed binary.

The Python package now contains many source-matched components and bounded
compositions, including membership updates, relation refresh, growth/thinning,
atlas-prior sampling/registration, and selected continuing-iteration slices.
They are not yet a complete end-to-end LesionTOADS workflow: unported wrapper
paths and options, the complete iteration/stop sequence, and full-output
validation remain open. The current helpers therefore provide no full
Lesion-TOADS segmentation-equivalence claim. Some helpers port dormant Java
methods or source branches that the selected wrapper does not call. See
[`lesion-toads-port-plan.md`](../lesion-toads-port-plan.md) for the remaining
source-porting and synthetic-validation steps.
