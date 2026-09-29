# AlgorithmTopologyCorrection Java reference phantom

These committed synthetic cases run the selected VandyCRUISE
`AlgorithmTopologyCorrection` source in a headless harness. All use scalar-image
mode, 3-D 6/18 connectivity, levels 0 and 1, min/max normalization, minimum
distance 0.05, and downward `object->background` propagation.

| Fixture | Shape | Coverage |
|---|---:|---|
| Root fixture (`spatial_hash`) | 9×9×9 | Many distinct intensity ranks, zero corner, one 1.0 seed and adjacent 0.5 paint value to exercise the strict `value > 0.5` adapter threshold. |
| `plateau/` | 9×9×9 | Broad repeated levels plus input values below 0 and above 1, testing plateau ties and normalization/clamping. |
| `boundary/` | 9×9×9 | Painted seed touches the image boundary. |
| `narrow_branch/` | 11×9×9 | Thin multi-voxel painted branch and adjacent threshold-edge paint value. |
| `disconnected_seeds/` | 9×9×9 | Two disconnected painted seeds, including a threshold-edge neighbor. |

Each `java_reference.nii.gz` is written by the original Java algorithm. Python
matches all five output arrays exactly: 3,807 voxels total, zero differing
values. Trace mode also matches all queue-pop coordinates and float32
priorities, including the number of topology-critical rejected candidates.
Across the five cases, event counts are 1,085, 909, 1,059, 373, and 1,222.
These are small harness comparisons, not a claim of whole-pipeline or clinical
equivalence.

Source copies, minimal headless MIPAV shims, and the harness are isolated under
`tests/java_reference/topology_correction/`. The source-copy hashes are:

| Java source | Selected checkout SHA-256 | Harness copy SHA-256 |
| --- | --- | --- |
| `AlgorithmTopologyCorrection.java` | `9dd7cd0db8c83e47a21059fe5acc52c7737cb0001222d6019619bf62526159a8` | `9dd7cd0db8c83e47a21059fe5acc52c7737cb0001222d6019619bf62526159a8` |
| `TopologyPropagation.java` | `781755f91a244569b0d3bb025aa201752996150f7d674ae7efc8de931d9475c9` | `c2c98157035a9c077e5de3a242cecadb9ce7106dd9df472d7a2ca5ba1e611c78` |
| `NeighborConnectivity.java` | `1021a8e76a0d57be16b473445e70c7c7b35cfeea8fa7805a78e2870f11c8e2db` | `1021a8e76a0d57be16b473445e70c7c7b35cfeea8fa7805a78e2870f11c8e2db` |
| `BinaryTree.java` | `2aa084673fc328702c501364edc27a8e38632c28c238ee332dfd5443df43d21b` | `2aa084673fc328702c501364edc27a8e38632c28c238ee332dfd5443df43d21b` |

The `TopologyPropagation.java` harness copy differs only by conditional
diagnostic tracing and a local boolean for the existing topology decision;
trace output is disabled in ordinary fixture generation. `ObjectProcessing`
signatures are a compile-time stub for the unreachable `inputSkel="intensity"`
branches. The exercised `paint_mask` path calls the copied
`AlgorithmTopologyCorrection`, `TopologyPropagation`, `NeighborConnectivity`,
and `BinaryTree` implementations.

Regenerate all fixture NIfTIs from the repository root with
`scripts/generate_java_topology_correction_reference.sh`, then report the
per-case Python comparison with
`PYTHONPATH=src .venv/bin/python tests/java_reference/topology_correction/compare.py`.
Trace a particular case by setting `TOPOLOGY_TRACE=1` and
`TOPOLOGY_TRACE_CASE=spatial_hash|plateau|boundary|narrow_branch|disconnected_seeds`
on the generation command, then pass the same case name to
`PYTHONPATH=src .venv/bin/python tests/java_reference/topology_correction/trace_compare.py`
to summarize queue events; add `--verbose` to inspect initialization and the
first queue decisions. Trace files are generated under ignored
`build/java-reference/`.

The Java `ModelImage` linear index has x as its fastest axis. The converter
reshapes flat Java values in Fortran order so the committed NIfTI coordinate
axes preserve the original x/y/z locations.

The same generator also writes the Dura two-call composition fixture under
`tests/data/synthetic/dura_topology_two_pass/`; its source and limits are
documented in that directory's README.

The shims provide only the in-memory image and progress APIs used by this run;
they do not reproduce MIPAV display, metadata, or file-I/O behavior.
