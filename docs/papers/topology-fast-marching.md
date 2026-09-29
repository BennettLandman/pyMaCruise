# Fast-marching topology correction: paper-to-code map

## Published method

Pierre-Louis Bazin and Dzung L. Pham, “Topology correction of segmented
medical images using a fast marching algorithm,” *Computer Methods and Programs
in Biomedicine*, 88(2), 182–190 (2007).
[DOI: 10.1016/j.cmpb.2007.08.006](https://doi.org/10.1016/j.cmpb.2007.08.006) ·
[PubMed record](https://pubmed.ncbi.nlm.nih.gov/17942182/) ·
[full text (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2128043/).

The paper describes topology propagation on a scalar or fuzzy segmentation,
using critical-point analysis and a fast-marching process so the requested
topology is respected across the field's level sets. This is different from
changing only one binary surface. This document records the conceptual mapping
and implementation boundaries; it does not claim that this Python module is a
standalone reimplementation of every option or experiment in the paper.

## Code mapping

| Published concept | Python implementation | Mapping and boundary |
| --- | --- | --- |
| Scalar/implicit segmentation field | `topology_correct` | The input volume is copied to float32, its level range is normalized according to the selected VandyCRUISE adapter behavior, and the field is propagated. The input remains a scalar field rather than being reduced to one binary surface. |
| A seed object from which propagation starts | `topology_correct` paint-mask initialization and `_downward_topology_propagation` | Paint voxels seed the object at the high level. This is the adapter's `paint_mask` / `object->background` workflow. The paper's method also permits other setup choices; the API here implements the selected wrapper path. |
| Fast-marching / topology-propagation order | `_downward_topology_propagation` | A priority heap processes the current descending front by scalar value. The source's critical/frontier/background/object states are represented by the module label constants. Candidate values and six-face neighbors schedule/re-schedule the front. |
| Reject topologically critical additions | `_is_regular_object_addition` and `topology_numbers_618` | An attempted addition is accepted only when the local foreground and background component counts are each one. A failed test marks the candidate `_CRITICAL`, preserving the topology of the propagated object. |
| Desired topology through digital connectivity | `_FOREGROUND_6`, `_BACKGROUND_18`, `topology_numbers_618` | This port uses the VandyCRUISE adapter's 3-D foreground-6/background-18 setting and excludes the neighborhood center from component counts, matching its `NeighborConnectivity` implementation. The paper's general method is not restricted to this one connectivity choice. |
| Field bounds and regularization | `topology_correct` arguments and `_downward_topology_propagation(min_distance=...)` | `lowest_level=0`, `highest_level=1`, `regularize=1e-4`, and `paint_threshold=0.5` are the selected wrapper defaults. They are API/source-adapter choices, not universal parameter values claimed by the publication. `regularize` becomes a value-space decrement scaled by the input range. |
| Domain boundary and propagation stop | `_downward_topology_propagation` | A one-voxel zero frame and stopping at the selected lower level reproduce the source adapter's behavior. They are implementation details to validate against the Java source, not paper-level invariants. |

The optional `multi_atlas_wm` branch in `topology_correct` performs a separate
V3 wrapper refinement after propagation: one six-neighbor minimum-filter pass,
then promotion of low corrected values in the refined atlas WM region. That
post-processing is not attributed to Bazin & Pham's fast-marching method.

### Source-code anchors used for this mapping

The selected VandyCRUISE Java source path is
`VandyCRUISE/src/edu/jhmi/rad/medic/libraries/TopologyPropagation.java`:
`initDownwardPaintLabels` initializes the paint-seeded state,
`propagateDownwardExactSmoothing` runs the downward queue propagation, and
`isRegularObjectPoint` applies the local regularity check. The connectivity
branches call `NeighborConnectivity.positiveEqual618Neighborhood3D` and
`negative18Neighborhood3D`. The adapter is
`VandyCRUISE/src/edu/jhu/ece/iacl/plugins/topology/MedicAlgorithmTopologyCorrection.java`,
which supplies the threshold, connectivity pair, levels, and regularization
parameters. The Python names above map to these behaviors, while selected
boundary details and argument defaults remain implementation-specific.

## Distinction from Han et al.'s GTCA

Han et al.'s “Topology correction in brain cortex segmentation using a
multiscale, graph-based algorithm” (*IEEE Transactions on Medical Imaging*,
21(2), 109–121, 2002) describes the graph-based topology correction algorithm
(GTCA). See the [PubMed record](https://pubmed.ncbi.nlm.nih.gov/11929099/) and
the authors' [GTCA project page](https://iacl.ece.jhu.edu/static/gtca/).

GTCA is a separate algorithm: its graph-based, multiscale optimization and cut
operations are not implemented in `topology.py`. The fact that Bazin & Pham
compare their propagation method with GTCA does not make the two methods
interchangeable. The local connectivity test here is part of the fast-marching
method and should not be described as GTCA.

## Scope of this mapping

This map explains why the code is structured as a scalar-field propagation
with local topology tests. Exact reproduction claims remain bounded by the
selected VandyCRUISE Java implementation and wrapper settings, and must be
supported by the project's Java-harness validations. In particular, the
published method's generality must not be used to imply that this adapter
supports every propagation direction, dimension, connectivity, or image type
available in the source algorithm.
