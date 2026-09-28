# Synthetic coarse-scale registration kernel images

These committed float32 probability-field NIfTIs exercise the selected TOADS
rigid-registration energy and coefficient formulas at pyramid scales 2 and 4.
They contain no clinical or subject data.

| File | Shape | RAS spacing | Definition |
|---|---:|---:|---|
| `memberships.nii.gz` | 13×12×11×3 | 0.8×1.2×1.5 mm | `((7*x + 5*y + 3*z + 2*k) mod 17) / 16` |
| `atlas_shapes.nii.gz` | 17×16×15×3 | 0.9×1.1×1.3 mm | `((2*x + 5*y + 3*z + 4*k) mod 11) / 10` for `x < 8`, otherwise 0. |

The fourth axis is class, with class 0 excluded by the registration objective.
Tests check positive RAS orientation and exact synthetic values. The committed
`java_reference.json` was produced by
`tests/java_reference/ToadsScaledRegistrationKernelsReferenceHarness.java`, a
standalone formula excerpt of the selected source mapping, interpolation,
derivative, and registration accumulation methods. The harness block-averages
each full-grid input in Java float32 loop order before evaluating scales 2 and
4. It records positive-prior and skipped-zero counts, which the Python test
checks against its scalar sampler. This is source-formula reference evidence,
not execution of the complete MIPAV/JIST atlas class.
