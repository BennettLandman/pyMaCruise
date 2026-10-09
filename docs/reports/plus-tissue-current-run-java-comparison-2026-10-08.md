# Held-fixed plus tissue-fusion Java source comparison — 2026-10-08

## Scope

This compares the preserved selected `CombineMAandTOADS.java` source class directly against Python `combine_multi_atlas_toads` on the six exact crop-grid inputs from the saved folded synthetic plus run. The grid is 66×76×66 and RAS. Java executes the original selected class through headless JIST parameter stubs, including its `MASImorph.Minimize3DFilter` call, default `ErodeMethod=3D`, and `GMconfident=0.001`.

This held-fixed comparison starts after upstream TOADS/LesionTOADS and registration. It does not claim whole-plus workflow equivalence or real-MRI validation.

## Reproduction

```sh
./scripts/compare_plus_tissue_update_saved_java.sh
```

The command records hashes for the workflow manifest and all six crop inputs in `build/plus-tgdm-current-java-2026-10-08/topology-provenance.json`, then compiles and runs the selected Java source class with small headless JIST stubs. The output order is converted back from the Java array traversal before comparing voxel values.

## Result

All three outputs match exactly across all **331,056 voxels**:

| Output | Changed voxels | Maximum absolute error | Mean absolute error |
|---|---:|---:|---:|
| Fused GM | 0 | 0 | 0 |
| Fused WM | 0 | 0 | 0 |
| Fused mask | 0 | 0 | 0 |

This is exact selected-source evidence for plus tissue fusion on a workflow-derived folded synthetic case. The standalone committed 7³ fixture continues to cover threshold equality and additional branch conditions; this larger comparison checks the composed real workflow inputs.

## Interpretation

The MaSeg/TOADS tissue-fusion adapter matches the selected Java class on the current synthetic plus run. Combined with the separate exact topology, MaSeg-WM erosion/refinement, and TGDM comparisons, this narrows the source-comparison gap through the plus surface stage. Registration, LesionTOADS source parity, complete plus workflow equivalence, and human-MRI validation remain open.
