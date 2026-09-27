# GVF port progress

The selected workflow is `run_GVF_jist.m` →
`MedicAlgorithmGradVecFlow` → VandyCRUISE `GradVecFlowGeneralized`. The
MaCRUISE wrapper supplies a scalar ACE enhanced-GM image, `Generalized`, lambda
1.0 by default (or configured `GVF_Lambda`), 200 iterations, and `normalize=false`.

The generalized implementation min/max normalizes its scalar input, computes
the source central gradient (forward/backward differences on volume edges),
forms `g = exp(-|grad f|²/lambda²)/6` and `h = 1-g`, and solves three
component-wise Poisson systems with the same coefficient fields. The source
uses a custom multigrid V-cycle and red/black Gauss-Seidel smoother. The Python
implementation will use a matrix-free SciPy sparse solver for that same
discretized equation; its convergence and component-wise residual must be
compared against Java before P28 can be validated.

NIfTI outputs will be RAS-oriented `.nii.gz` three-component vector fields;
non-RAS input must be normalized with a warning. The 3-vector image
representation uses shape `(x, y, z, 3)` with components in index-axis order,
and physical-coordinate component convention will be confirmed against the
source `ImageData`/TGDM use before downstream P29 integration.

Synthetic source-reference fixtures and solver checks have not yet been
implemented. No actual imaging scans are part of this stage.
