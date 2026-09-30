# P18czm synthetic registration handoff

This P18czm handoff fixture deliberately reuses the committed synthetic
P18czl input volumes under `../toads_lesion_atlas_type_membership/` rather than
duplicating them. The subject is a folded, asymmetric 9×8×7 two-channel RAS
phantom; its atlas shape maps occupy a separate 11×10×9 template grid, with 13
classes. Both grids use `.nii.gz` NIfTI files. The membership maps supplied to
registration are taken from the pinned Java output in that fixture, including
LESION at class index 10.

`java_reference.json` was emitted by
`tests/java_reference/ToadsLesionAtlasRegistrationHandoffReferenceHarness.java`
from the bundled lesion-TOADS descriptor. It records the selected wrapper call
order, source identifiers, and the descriptor-ordered class indices passed to
registration. Numeric registration kernels have separate selected-source
Java formula harnesses documented in `docs/toads-atlas-registration.md` and
`docs/reports/toads-rigid-registration-update-2026-09-29.md`; the focused
Python test feeds the Java-generated P18czl membership bits into the
registration handoff and verifies that its full transform/trial state equals
direct invocation of the already-compared registration API.

All images are deterministic synthetic data. No human imaging data or
external source checkout is required to run the P18czm tests.
