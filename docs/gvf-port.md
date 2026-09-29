# GVF port progress

The selected workflow is MaCRUISE `run_GVF_jist.m` → JIST
`MedicAlgorithmGradVecFlow` → VandyCRUISE `GradVecFlowGeneralized`. The base
wrapper defaults to `Generalized`, lambda 1.0, 200 maximum iterations, and
`normalize=false`. MaCRUISEplus overrides lambda to 0.2 while retaining 200
iterations and no output normalization. `run_plus_gvf_stage` in
`src/pymacruise/plus_pipeline.py` forwards those plus settings to the existing
Python GVF adapter; it does not implement another solver.

## Source behavior and Python implementation

The Java path min/max normalizes the scalar volume, computes central
differences with one-sided differences on the image boundary, and obtains the
gradient magnitude with JIST `ImageDataMath.vectorMagnitude4D` (square root of
the sum of squared components). Despite the Java local variable name
`SqrMagf`, that magnitude is used directly in `exp(-magnitude/lambda²)`;
`g=exp(...)/6` and `h=1-exp(...)`. It then solves three systems
`h*u - g*laplacian(u) = h*gradient_component` using a multigrid V-cycle with
red/black Gauss-Seidel smoothing. The earlier Python implementation squared
the magnitude in the exponential; Java reference comparison caught and fixed
that discrepancy.

`src/pymacruise/gvf.py` solves the same discrete mirrored-boundary system with
matrix-free SciPy conjugate gradients. `generalized_gvf_image` reorients input
images to RAS with a warning when needed and returns a component-last float32
vector NIfTI with RAS affine/header. Components 0, 1, and 2 follow array axes
x, y, and z. This matches VandyCRUISE `Gradient3d` and the downstream TGDM
advection code, which pairs those components with x/y/z finite differences.
The selected legacy equations are voxel-grid based and do not apply physical
voxel-spacing factors; this preserves source behavior. Synthetic fixtures use
1 mm isotropic grids.

The constant-image edge case is handled deliberately: Java's `(f-min)/(max-min)`
produces NaNs when the image is constant, while Python returns a finite zero
field. This defensive behavior is covered and documented as a source edge-case
deviation.

## Synthetic Java comparison

The test-only harness compiles the selected VandyCRUISE GVF sources against
minimal JIST/framework stubs; no MIPAV or JIST application is required. Source
files retain their author comments and CVS revision metadata. Reproduction
instructions and fixture provenance are in
[`tests/java_reference/README.md`](../tests/java_reference/README.md) and
[`tests/data/synthetic/gvf/java_reference/README.md`](../tests/data/synthetic/gvf/java_reference/README.md).

The committed 16×12×10 RAS `.nii.gz` fields include a mixed-axis ramp, an
asymmetric step, and overlapping asymmetric Gaussian blobs. The Java
references use lambda 1.0, up to 200 iterations, and no vector normalization.
The generated `spec.json` reports:

| Synthetic input | Mean absolute component error | Maximum absolute component error | RMSE |
|---|---:|---:|---:|
| Mixed-axis ramp | 1.30e-8 | 9.69e-8 | 1.71e-8 |
| Asymmetric step | 6.69e-6 | 3.26e-5 | 1.04e-5 |
| Asymmetric Gaussian blobs | 1.17e-6 | 8.55e-6 | 1.80e-6 |

Acceptance for these fields is mean absolute error ≤1e-5 and maximum absolute
error ≤5e-5. Analytic gradient, solver residual, lambda response, finite-value,
RAS vector intent, and non-RAS reorientation tests also pass. The latter
compares a physically equivalent flipped input to its RAS counterpart and
confirms component signs follow the reoriented axes. No actual imaging data is
used in this step.

Regenerate golden outputs and their checksums with:

```sh
python scripts/generate_java_gvf_reference.py
```

This requires a JDK and the project Python dependencies. Temporary classes
and binary interchange files are kept under ignored `build/`.
