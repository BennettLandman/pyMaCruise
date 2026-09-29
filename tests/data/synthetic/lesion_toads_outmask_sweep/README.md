# LesionTOADS OUTMASK volume fixture

The `5×4×6` synthetic float32 volume has two ordinary OBJECT classes and one
OUTMASK class, per-voxel class/outmask data terms, an outlier membership map
with zero, interior, and unit values, mixed direct/prior relations, and one
zero-shape OUTMASK voxel. The selected-Java formula harness follows source
x/y/z in-place order, applies six-face smoothing, blends class and outmask
memberships, updates the outlier map, and normalizes at each voxel. Committed
references include all interior class and outlier membership bits.

This fixture validates the OUTMASK smoothing/blending/normalization branch
only; atlas priors and mixed lesion/black-hole maps are not included.
Regenerate with `scripts/generate_java_lesion_toads_outmask_sweep_reference.sh`.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
