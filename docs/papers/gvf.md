# Generalized gradient vector flow (GVF)

## Published model

Xu, C. and Prince, J. L. “Generalized gradient vector flow external forces for
active contours.” *Signal Processing*, 71(2), 131–139, 1998.
[DOI: 10.1016/S0165-1684(98)00140-6](https://doi.org/10.1016/S0165-1684(98)00140-6).
[Author-hosted paper PDF](https://iacl.ece.jhu.edu/pubs/p105j.pdf).

The paper introduces GGVF as a vector diffusion model. At equilibrium, each
spatial component of the force field \(v\) satisfies the same equation:

\[
v_t = g(|\nabla f|)\nabla^2 v - h(|\nabla f|)(v-\nabla f),
\qquad v_t=0 \text{ at equilibrium}.
\]

Here \(f\) is an edge map, the first term smooths the field, and the second
pulls it toward the edge-map gradient. In the paper's example weighting pair,
the smoothing weight is an exponential of the **squared** edge-map-gradient
magnitude and the data weight is one minus that exponential (Eqs. 8–9). These
model-level definitions describe the method; they do not prescribe the
VandyCRUISE wrapper's parameter values or its exact discrete solver.

Using \(s=|\nabla f|\), that example pair is
\(g_{\mathrm{paper}}(s)=\exp(-s^2/K)\) and
\(h_{\mathrm{paper}}(s)=1-g_{\mathrm{paper}}(s)\), where \(K\) controls the
tradeoff between smoothing and gradient conformity.

## Mapping to this implementation

| Paper or source behavior | Python mapping |
| --- | --- |
| The input edge map is differentiated to provide the data vector field. | `generalized_gvf` min/max-normalizes the scalar volume, then _source_gradient applies central differences internally and one-sided differences at the volume boundary. This follows VandyCRUISE `GradVecFlowGeneralized.solve` and `Gradient3d`; it is a voxel-grid stencil with no physical-spacing factors. |
| GGVF applies shared, spatially varying weights to three independently diffused vector components. | The magnitude-derived `edge_weight`, `g`, and `h` arrays are shared by the three component solves. Each component's right-hand side is `h * gradient_component`. |
| The equilibrium PDE balances smoothing against attraction to the input gradient. | The matrix-free operator and right-hand side encode `h*u - g*Delta_6(u) = h*gradient_component`, rearranged from the selected source equilibrium equation. _source_laplacian is the six-neighbor stencil and clamps boundary indices to the nearest in-volume voxel, matching `Poisson3d.residue`. |
| The paper describes a model, not this package's array or file convention. | `generalized_gvf` returns the three components on array axes 0, 1, 2. `generalized_gvf_image` first reorients scalar NIfTI input to RAS (warning when required) and writes a component-last vector NIfTI with that RAS affine. |

### VandyCRUISE-specific parameterization and numerical method

The selected VandyCRUISE source is GradVecFlowGeneralized with its helper
Poisson3d, called by MedicAlgorithmGradVecFlow. It computes the ordinary
Euclidean magnitude \(m=|\nabla f|\) using vectorMagnitude4D (even though the
local Java variable is named SqrMagf), then sets

\[
q=\exp(-m/\lambda^2),\qquad
g_{\mathrm{source}}=q/6,\qquad
h_{\mathrm{source}}=1-q.
\]

Thus this selected implementation is **not a literal transcription** of the
paper's example weighting formula: it uses \(m\), rather than \(m^2\), in the
exponent, and divides the smoothing coefficient by six for the source's
six-neighbor operator. Those are deliberate source-compatibility details,
not claims about the general GGVF model. The active MaCRUISE wrapper selects
the Generalized method, sets lambda=1 and a 200-iteration maximum, and
requests unnormalized vectors; alternative pipeline profiles may set another
lambda (for example 0.2).

VandyCRUISE's Poisson3d solves the component equations with a multigrid
V-cycle and red/black Gauss-Seidel smoothing. Python applies the same selected
discrete system matrix-free but solves each component with SciPy conjugate
gradients (scipy.sparse.linalg.cg). This is a solver substitution: the
operator, input weighting, and boundary stencil are carried over, but solver
trajectories and finite-precision results are not identical. Python applies
the operator in float64 and stores the returned field as float32; its
iteration cap and relative tolerance are explicit function arguments.
Agreement claims therefore refer to the recorded synthetic Java harness
comparison, not to identity of algorithms or exact arithmetic.

## Code pointers

- [gvf.py](../../src/pymacruise/gvf.py): input normalization, source weights,
  per-component system, SciPy CG substitution, RAS vector-image wrapper, and
  finite-difference stencils.
- VandyCRUISE reference implementation: GradVecFlowGeneralized.java,
  Gradient3d.java, Poisson3d.java, and the JIST adapter
  MedicAlgorithmGradVecFlow.java in the supplied VandyCRUISE source tree.
- [gvf-port.md](../gvf-port.md): source behavior, synthetic Java harness,
  measured errors, and validation scope.

The Xu–Prince paper formulates GGVF for 2-D active contours. The VandyCRUISE
code adapts the model to 3-D medical-image volumes; this implementation follows
that 3-D source path rather than claiming the paper itself specifies the 3-D
stencil or source parameterization.
