# Generalized GVF Java reference fixtures

These synthetic 1 mm isotropic fields are stored as RAS `.nii.gz` images.
`*_input.nii.gz` files are analytic scalar inputs; `*_gvf_java.nii.gz` files
are outputs from the copied VandyCRUISE `GradVecFlowGeneralized` Java source
using lambda 1.0, 200 maximum iterations, and vector normalization disabled.
Cases include a mixed-axis ramp, an asymmetric thresholded step, and two
overlapping asymmetric Gaussian blobs.

The fixture generator compiles the supplied Java algorithm and runs it using
the minimal test-only JIST stubs; it also records componentwise Python-vs-Java
MAE, maximum absolute error, and RMSE in `spec.json`.

```sh
python scripts/generate_java_gvf_reference.py
```

A JDK and the project Python dependencies are required. Transient Java class
files and binary interchange files are written under ignored `build/`.
Source provenance and the vector component convention are summarized in
`tests/java_reference/README.md` and `docs/gvf-port.md`.
