# AlgorithmTopologyCorrection Java reference phantom

This deterministic 9×9×9 synthetic case runs the selected VandyCRUISE
`AlgorithmTopologyCorrection` source in a headless harness. The source image is
a deterministic spatial-hash intensity field with a zero corner. The paint
image has one 1.0 seed and an adjacent 0.5 voxel to exercise the JIST adapter's
strict `value > 0.5` rule. It uses scalar-image mode, 3-D 6/18 connectivity,
levels 0 and 1, min/max normalization, minimum distance 0.05, and downward
`object->background` propagation.

`java_reference.nii.gz` is the output from the original Java algorithm. The
current Python output matches it exactly on this fixture (0 of 729 differing
voxels; maximum absolute difference 0). In trace mode, all 1,085 popped queue
coordinates and float32 priorities also match. This is a small harness
comparison, not a claim of whole-pipeline equivalence.

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

Regenerate all three NIfTI files from the repository root with
`scripts/generate_java_topology_correction_reference.sh`, then report the
Python comparison with
`PYTHONPATH=src .venv/bin/python tests/java_reference/topology_correction/compare.py`.
Add `TOPOLOGY_TRACE=1` to the generation command, then run
`PYTHONPATH=src .venv/bin/python tests/java_reference/topology_correction/trace_compare.py`
to inspect initialization and queue decisions. Trace files are generated
under ignored `build/java-reference/`.

The shims provide only the in-memory image and progress APIs used by this run;
they do not reproduce MIPAV display, metadata, or file-I/O behavior.
