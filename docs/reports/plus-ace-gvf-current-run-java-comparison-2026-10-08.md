# Held-fixed plus ACE and GVF Java comparison — 2026-10-08

## Scope

This compares selected MaCRUISEplus ACE mode 3 and VandyCRUISE generalized GVF against Python on reconstructed inputs from the saved 64×80×64 folded synthetic plus run. The crop-ready grid is 66×76×66 in RAS. The Python replay first checks all six saved TGDM-stage products exactly against the workflow outputs; ACE/GVF inputs are then generated using the same plus options and crop inputs.

## Reproduction

```sh
./scripts/compare_plus_ace_gvf_saved_java.sh
```

The Java ACE run uses the source after-Gaussian CSF merge (`place=3`), `alpha=1.8`, and the selected threshold/diffusion settings. GVF uses `lambda=0.2`, 200 iterations, and no vector normalization. Both run under an 8 GB RSS watchdog. The sampled peaks were **0.265 GB for ACE** and **0.249 GB for GVF**.

## Results

All three ACE outputs are voxel-exact across 331,056 voxels: enhanced GM, cleaned skeleton, and thinned skeleton have zero changed voxels and zero absolute error.

GVF differs only at floating-point rounding scale: mean absolute error **2.30×10⁻⁸**, maximum absolute error **2.13×10⁻⁷** across 993,168 vector components. This passes the established Java-reference limits (mean ≤1×10⁻⁵, maximum ≤5×10⁻⁵).

## Interpretation

This adds selected-source evidence for the plus ACE and GVF stages on the current workflow-derived synthetic case. Along with exact tissue fusion, topology correction/refinement, TGDM, and JIST surface extraction, it validates the current synthetic plus surface chain through its major algorithm stages with upstream registration and LesionTOADS held fixed. Full composite parity and human-MRI validation remain open.
