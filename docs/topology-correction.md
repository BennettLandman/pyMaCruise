# Topology Correction

**Status:** The selected Java-class output now matches five committed synthetic fixtures (3,807 voxels total) exactly, and trace mode matches 4,648 queue pops including critical-candidate counts. This validates the selected low-level paint-mask algorithm cases; it does not validate whole-workflow or clinical-data equivalence.

## Selected MaCRUISE behavior

The corrected V3 workflow invokes the JIST `MedicAlgorithmTopologyCorrection` adapter with a scalar WM image, paint mask, `(6,18)` foreground/background connectivity, `object->background` propagation, intensity initialization, normalized levels 0–1, min/max normalization enabled, and regularization `1e-4`. The CLI wrapper has no stop parameter. In `AlgorithmTopologyCorrection`, the field is padded by one zero voxel on every face before deriving its min/max; the paint mask uses the adapter's `>0.5` threshold. The algorithm sets paint voxels to the high level and runs downward exact smoothing; that propagation method enables threshold stopping at the low level.

`run_TopoCorrect.m` then performs the selected V3 multi-atlas refinement. It erodes the atlas WM with the source's default 3D center-plus-six-neighbor integer minimum filter (the output border remains zero), then changes corrected WM values `<200` to `255` wherever the eroded atlas WM is positive. The optional Python `multi_atlas_wm` argument accepts the original, un-eroded atlas volume and performs both operations.

## Python implementation

`src/pymacruise/topology.py` implements the default scalar `(6,18)` paint-mask pathway using NumPy arrays, SciPy connected-component labeling, and a max-priority heap. A candidate joins the object only when its local foreground 6-component and background 18-component counts are each one. The returned corrected field is float32 and retains the source image's RAS grid. The atlas-refined field is a separate result, like the two products in the MATLAB wrapper.

The connectivity counts follow the supplied VandyCRUISE `NeighborConnectivity` implementation: foreground connectivity is six-connected with corner voxels excluded; background connectivity uses the 18-neighborhood. The propagation and local connectivity checks are implemented directly rather than replaced with binary closing, hole filling, or a generic topology filter.

## Synthetic validation

`tests/test_topology.py` checks:

- Exact local connectivity numbers for regular and critical neighborhoods.
- Two distinct painted components remain disconnected under downward propagation.
- A seeded torus retains cubical Euler characteristic 0 (one handle).
- A seeded hollow shell retains cubical Euler characteristic 2 (one enclosed cavity).
- The V3 atlas erosion/refinement runs as a separate output and preserves RAS grid metadata.
- Non-RAS image/paint inputs warn and reorient; empty masks, non-finite image values, and mismatched grids fail clearly.
- Repeated runs are deterministic on the synthetic cases.

The first headless Java run exposed an initial implementation mismatch of 621/729 voxels. Source-level follow-up corrected face-connected foreground counting, strict tie ordering in VandyCRUISE's custom `BinaryTree`, threshold-stop fill behavior, and float32 operation order. The five-case comparison now matches 3,807/3,807 voxels exactly. Cases exercise spatial-hash distinct ranks, broad plateaus and input values outside 0–1, boundary-touching paint, a narrow painted branch, and disconnected seeds; trace mode also matches all 4,648 queue events and critical rejection counts exactly. The committed fixtures, isolated selected-source harness, regeneration command, hashes, and limits are recorded in [the P25 comparison report](reports/topology-correction-java-comparison-2026-09-29.md). The V3 atlas refinement is separately tested; neither result establishes full-workflow or clinical-data equivalence.
