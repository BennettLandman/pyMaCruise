# P25 topology correction Java comparison

## Result

The committed 9×9×9 Dura-settings phantom now matches the selected
VandyCRUISE Java `AlgorithmTopologyCorrection` output exactly:

| Measure | Result |
|---|---:|
| Voxels compared | 729 |
| Different output voxels | 0 |
| Maximum absolute difference | 0 |
| Java / Python float64 sum | 126.464972 / 126.464972 |
| Java / Python queue pops | 1,085 / 1,085 |
| Queue coordinates and float32 priorities | Exact for all 1,085 pops |

The phantom uses scalar-image mode, 3-D `(6,18)` connectivity, downward
`object->background` propagation, normalized levels 0–1, min/max normalization,
paint threshold `>0.5`, and regularization `0.05`. It contains a deterministic
spatial-hash intensity field with a zero corner, one 1.0 seed, and one adjacent
0.5 paint voxel. Inputs and the Java result are under
`tests/data/synthetic/topology_correction/`.

## Source and debugging trail

The isolated harness compiles selected copied source files for
`AlgorithmTopologyCorrection`, `TopologyPropagation`, `NeighborConnectivity`,
and `BinaryTree` with small in-memory MIPAV/UI shims. `ObjectProcessing` is a
compile-time stub for an unreachable intensity-skeleton branch; the tested
paint-mask branch executes the copied algorithm sources. The Java output is
written by that run and converted to NIfTI only for voxelwise comparison.

The first Python implementation comparison differed at 621/729 voxels
(maximum absolute difference 0.411793947). Tracing the source exposed four
implementation details required for this fixture:

1. The foreground component count includes only 6-connected components that
   touch one of the candidate's six face neighbors; edge-only components do
   not contribute.
2. Equal-priority queue operations follow VandyCRUISE's custom
   `BinaryTree` strict-comparison insertion/removal order, rather than the
   tie ordering of Python's standard heap.
3. With the adapter's threshold stop enabled, unprocessed background and
   critical labels are filled with the configured low level, not the final
   popped queue priority.
4. Normalization, minimum-distance updates, and queue priorities preserve the
   source's float32 operation points.

After applying these source-derived rules, all queue events and final values
match on this fixture. The initial mismatch is retained here as debugging
history, not as the current result.

## Reproduction

From the repository root:

```sh
scripts/generate_java_topology_correction_reference.sh
PYTHONPATH=src .venv/bin/python tests/java_reference/topology_correction/compare.py
```

The first command needs a JDK (`javac`/`java`) and regenerates the small
reference fixtures. The second prints the five image metrics above. Selected
checkout and isolated harness-copy SHA-256 values:

| File | SHA-256 |
|---|---|
| `AlgorithmTopologyCorrection.java` | Selected and harness: `9dd7cd0db8c83e47a21059fe5acc52c7737cb0001222d6019619bf62526159a8` |
| `TopologyPropagation.java` | Selected: `781755f91a244569b0d3bb025aa201752996150f7d674ae7efc8de931d9475c9`; harness (diagnostic trace instrumentation): `c2c98157035a9c077e5de3a242cecadb9ce7106dd9df472d7a2ca5ba1e611c78` |
| `NeighborConnectivity.java` | Selected and harness: `1021a8e76a0d57be16b473445e70c7c7b35cfeea8fa7805a78e2870f11c8e2db` |
| `BinaryTree.java` | Selected and harness: `2aa084673fc328702c501364edc27a8e38632c28c238ee332dfd5443df43d21b` |

The `TopologyPropagation.java` copy adds conditional diagnostic tracing for
queue events and decisions; ordinary generation disables trace mode. The
instrumentation and provenance are described in the fixture README.

## Scope and remaining P25 cases

This is a selected Java-class harness match on one compact phantom, not full
MaCRUISE, Dura workflow, or clinical equivalence. P25 remains Started. Add
independent synthetic fixtures for boundary-touching support, masked regions,
multiple distinct intensity levels, broad plateaus/tie priorities, narrow
branches, and critical-point decisions; compare both image outputs and queue
events where applicable. Then run and record the full pytest suite on the
committed implementation before marking any broader source-matching row Done.
