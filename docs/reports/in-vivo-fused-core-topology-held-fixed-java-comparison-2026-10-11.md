# Human fused-core WM-fill topology correction: held-fixed Java comparison

**Result: exact agreement on the saved human topology-corrected WM-fill stage.**
The selected Java source replay and Python output match exactly at every voxel
in the 141 × 173 × 131 RAS crop (3,195,483 voxels): zero unequal float32
values, maximum absolute difference 0, and mean absolute difference 0. The
five sampled superlevel masks are also identical (Dice 1.0 at each level).

This isolates the topology-correction operation. It does not make the human
base workflow an accepted result: upstream MaSeg tissue-support metrics remain
below the project’s 0.97 target, and this comparison does not validate the
full JIST workflow, refinement, surface extraction, or thickness.

## Inputs and settings

The comparison used these existing products from the corrected, fused human
base-core run:

- Scalar input: `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/wmfill_cropped.nii.gz`
- Paint mask: `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/wmmask_cropped.nii.gz`
- Saved Python result: `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/wmfill_topology_corrected.nii.gz`

All share the RAS grid, shape 141 × 173 × 131, and affine
`diag(1,1,1,1)` with translation `(-70.5, -85.5, -16.5)` mm. The strict
paint predicate is `>0.5`. Both implementations received scalar-image mode,
regularization `1e-4`, 3D, object/background connectivities 6/18, levels 0–1,
min/max scaling, and object-to-background propagation. This matches
`toads_cruise_core.py`.

Input SHA-256 values:

| File | SHA-256 |
| --- | --- |
| `wmfill_cropped.nii.gz` | `26e5e3baf0a73c4d5cfc0b9ed291e25bca6255c06e56a79ed19d4f6dd702e3aa` |
| `wmmask_cropped.nii.gz` | `df037cb35701a6d1bb18357b04ab05a3cfc14e171c48184d9464e337762a340a` |
| Python `wmfill_topology_corrected.nii.gz` | `477769a97df3f202f60d5d7c0555ed2a5768e59fdf44e5ae96d8b2b44374a5ef` |
| Replayed Java raw output (`build/.../topology-java-current-human/java-output.f32le`) | `ec3b5ba95e40711aa847709cd9f2b7140a46a8d65391a7b293e9a8f654060d83` |

## Comparison results

| Measure | Result |
| --- | ---: |
| Voxels | 3,195,483 |
| Unequal float32 values | 0 |
| Maximum / mean absolute scalar difference | 0 / 0 |
| Java / Python scalar range | `[0, 255.00003051757812]` / `[0, 255.00003051757812]` |
| Positive-field 6-connected components | 155 |

| Superlevel threshold | Foreground voxels, Java / Python | Differing voxels | Dice |
| ---: | ---: | ---: | ---: |
| 0 | 1,061,758 / 1,061,758 | 0 | 1.0 |
| 63.75 | 454,131 / 454,131 | 0 | 1.0 |
| 127.5 | 408,515 / 408,515 | 0 | 1.0 |
| 191.25 | 391,714 / 391,714 | 0 | 1.0 |
| 254 | 388,806 / 388,806 | 0 | 1.0 |

## Java harness and scope

The runner compiles and invokes the actual selected source files checked into
`tests/java_reference/topology_correction`:

- `upstream/edu/jhmi/rad/medic/algorithms/AlgorithmTopologyCorrection.java` —
  SHA-256 `9dd7cd0db8c83e47a21059fe5acc52c7737cb0001222d6019619bf62526159a8`
- `upstream/edu/jhmi/rad/medic/libraries/TopologyPropagation.java` —
  SHA-256 `c2c98157035a9c077e5de3a242cecadb9ce7106dd9df472d7a2ca5ba1e611c78`
- `upstream/edu/jhmi/rad/medic/libraries/NeighborConnectivity.java` —
  SHA-256 `1021a8e76a0d57be16b473445e70c7c7b35cfeea8fa7805a78e2870f11c8e2db`
- `upstream/edu/jhmi/rad/medic/structures/BinaryTree.java` —
  SHA-256 `2aa084673fc328702c501364edc27a8e38632c28c238ee332dfd5443df43d21b`

The `ModelImage`, `AlgorithmBase`, and UI classes are small in-memory framework
shims; the topology algorithm and its propagation/connectivity/heap code are
the selected source implementations. Therefore this does execute the actual
selected topology algorithm source, but it is **not** a full MIPAV application
or JIST workflow execution. The generated human harness and all raw data are
ignored build artifacts, not committed MRI samples.

## Reproduction and resource use

From the repository root:

```sh
cd /Users/landmaba/dev/mipav/pymacruise
PYTHONPATH=src .venv/bin/python scripts/compare_human_fused_topology_java.py
```

The runner compiles in the output directory, exports NIfTI arrays as
little-endian float32 in Fortran order (x-fastest, matching `ModelImage`), then
compares the Java raw output against the saved Python NIfTI. Java was OpenJDK
21.0.11, Python 3.13.2. Java heap was capped at 2,304 MiB; OS child-process
high-water RSS accounting reported 838,303,744 bytes (0.781 GiB), below the
3.8 GiB guard. Java execution took 1.91 s; compilation took 0.54 s. Machine
JSON and raw outputs are in the ignored `build/.../topology-java-current-human/`
directory; the runner writes exact hashes and metrics to `comparison.json`.
