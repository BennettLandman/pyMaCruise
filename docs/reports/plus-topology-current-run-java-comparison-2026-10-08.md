# Held-fixed plus topology correction Java comparison — 2026-10-08

## Scope

This compares selected VandyCRUISE `AlgorithmTopologyCorrection` against the Python topology correction on the crop-grid data derived from the saved 64×80×64 folded synthetic plus run. The inputs are 66×76×66 on an RAS grid. Python rebuilds the selected tissue update and routes its fused WM, brain paint mask, and raw MaSeg WM image into topology correction.

The Java harness compares the scalar topology correction, source `MASImorph.Minimize3DFilter` MaSeg-WM erosion, and the literal `Refine_Topo_by_MA.m` threshold/assignment on the held-fixed crop data. This does not establish complete plus-workflow equivalence or real-MRI validity.

## Reproduction

```sh
./scripts/compare_plus_topology_saved_java.sh
```

The generated input provenance records the exact saved crop-input hashes and workflow run signature. The harness runs selected Java `AlgorithmTopologyCorrection` with the wrapper's (6,18) connectivity, paint mask, scalar-image mode, min/max scaling, and `0.05` normalized minimum distance. The input arrays and output are compared in the source's voxel-storage order.

## Result

The Java run completed successfully. The corrected scalar field, eroded MaSeg WM mask, and final MaSeg-refined topology field all match Python **bit-for-bit across all 331,056 voxels**: zero changed voxels and zero maximum/mean absolute error. The corrected-field Dice at zero is 1.0. This is exact selected-source evidence on a larger workflow-derived folded synthetic input.

The RSS watchdog returned `0.000 GB` for this short JVM process. That value is below useful measurement resolution and is not reported as a peak. The command used an 8 GB RSS limit.

## Interpretation

This validates the selected topology-correction core and MaSeg-WM erosion/refinement on the saved plus inputs after Python tissue fusion. It does not validate the full tissue-fusion/topology composite against Java, upstream registration or LesionTOADS, or human MRI. Whole plus parity, upstream registration/LesionTOADS, and human-MRI validation remain open.
