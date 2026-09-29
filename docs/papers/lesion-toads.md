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
| [`lesion_atlas_priors.py`](../../src/pymacruise/lesion_atlas_priors.py) | Parses selected modality intensity rows and lesion/black-hole scalar priors from the text descriptor, retaining declared structure order. | Supplies numerical statistical-atlas and lesion-model inputs used by the paper's model. | It does not load topology or shape volumes, register an atlas, sample transformed priors, or construct segmentation memberships. Descriptor values and malformed-input behavior are covered by P39b's source comparison. |
| [`lesion_toads_preprocessing.py`](../../src/pymacruise/lesion_toads_preprocessing.py) | Clamps negative intensities, computes per-channel robust bounds, normalizes by robust span, and initializes lesion/black-hole centroids using the selected source formula. | Prepares multichannel intensity values and lesion-class initialization for the statistical segmentation framework. | Exact arithmetic/order is traced to selected VandyCRUISE `LesionToadPreprocess.java`; it is not a transcription of every preprocessing choice in the paper. P39a validates this helper against its Java harness on a synthetic fixture. Zero robust spans intentionally raise `ValueError` rather than reproducing Java NaN/infinity. |
| [`lesion_distance.py`](../../src/pymacruise/lesion_distance.py) | Computes the selected six-connected, capped voxel-step distance to a boolean mask. | A low-level distance support primitive associated with the selected implementation's propagation/topology machinery. | This helper is not the paper's topology-preserving fast-marching evolution, has no topology atlas or physical-spacing metric, and does not segment lesions. P39c validates only this primitive against a Java harness. |

## Version and completeness limits

The project selected VandyCRUISE as the source tree for the port, but the
exact deployed plugin version used by historical MaCRUISEplus runs is still
unknown. Consequently these mappings describe the cited 2010 method at a
conceptual level and the specific methods/data parsed from the selected
VandyCRUISE source; they do not establish that the selected tree is identical
to the paper implementation or the historical deployed binary.

The iterative lesion engine and full atlas path remain incomplete. In
particular, the Python package does not yet implement topology/shape-atlas
loading and transformed-prior sampling, lesion and black-hole memberships,
relation/boundary maps, topology-preserving competition/propagation, or final
Lesion-TOADS output assembly. The current helpers therefore provide no
Lesion-TOADS segmentation-equivalence claim. See
[`lesion-toads-port-plan.md`](../lesion-toads-port-plan.md) for the remaining
source-porting and synthetic-validation steps.
