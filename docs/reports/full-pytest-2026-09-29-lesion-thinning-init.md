# Full pytest evidence: LesionTOADS thinning frontier initialization

- **Implementation commit:** `7820fc2d22ecee9ea531fc3c8c05aa44303fe1fe`
- **Command:** `.venv/bin/pytest -q`
- **Result:** 600 passed, 2 expected RAS-reorientation warnings, 80.34 seconds.
- **Focused selection:** `tests/test_lesion_toads_thinning_initialization.py` — 2 passed.
- **Java reference:** selected `LesionToadSegmentation.java` excerpt, SHA-256
  `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`;
  the copied selected-source `BinaryHeap4D.java` SHA-256 is
  `9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a`.

The initializer matches both `useLesionMax` cases for class volumes,
per-class boundary counts, interior skeleton/availability state, every queue
insertion priority and payload, and the full equal-priority Java max-heap pop
sequence. A non-boundary interior voxel, untouched outer planes, lesion and
black-hole scoring, and the source `bestMembership()` running-maximum quirk
are exercised. Inputs and references are committed under
`tests/data/synthetic/lesion_toads_thinning_initialization/` as RAS `.nii.gz`
volumes and JSON numerical expectations.

This validates only `initThinningLabels()` and queue initialization; it does
not validate thinning queue events, topology-gated segmentation changes,
membership updates, or full LesionTOADS equivalence.
