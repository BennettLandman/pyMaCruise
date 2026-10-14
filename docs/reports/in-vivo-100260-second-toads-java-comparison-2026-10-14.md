# Case 100260 second-pass TOADS source comparison (2026-10-14)

## Result

The current Python second TOADS pass was compared with the selected Java
`ToadsWorkflowSourceHarness` using the exact same Dura-stripped T1, brain mask,
and source atlas. Both runs used the source-selected registration settings and
a 50-iteration limit. Their final classifications are **exactly equal**:
zero differing voxels across 10,813,440 and Dice 1.0 for each of the 11 labels.
Both outputs are on the shared 176×240×256 RAS grid.

Python reported iteration-limit termination at iteration 50 with energy ratio
`0.00145347673`, above the `1e-4` threshold. The Java harness also completed
iteration 50 and returned its final classification; its `completed=true`
status means the harness execution succeeded and does not certify algorithm
convergence. Thus this is strong source-equivalence evidence at the capped
checkpoint, while natural convergence remains an unresolved workflow blocker.

The Java classification is also identical to the earlier Python 100260
second-pass map. The current Python classification and its four Dura handoff
maps reproduce that prior run exactly. These results establish reproducible
selected-source behavior for this input, not a complete BASE Java workflow or
independent anatomical acceptance.

## Held-fixed input and provenance

- Dura-stripped T1: shape 176×240×256, RAS, spacing approximately
  1.00000155×1×1 mm, SHA-256
  `34ea8cb2999eb92f7c6638209951aae9cd2857013279be0ddcb6f8c7925c3382`.
- Same-grid Dura brain mask: SHA-256
  `ec12f6ce96fe21580e9a4eea472235905527fb7fa139139479249bedc4896b68`;
  1,349,513 positive voxels.
- The Java input serialization applies that mask to the T1, then writes
  little-endian float32 in MIPAV x-fastest order; its SHA-256 is
  `9702219547293bc4e4ce60174ae2ed39d54a2f1cc8899cfc42e2747205e19db5`.
- Atlas descriptor: the original VandyCRUISE
  `cruise-atlas-10obj-toads2008.txt`, not the converted Python bundle.
- Harness: `tests/java_reference/ToadsWorkflowSourceHarness.java`, SHA-256
  `397cb71fc3a8a18186fc132321bae54b7aa3691a9fd5331725f1b8287885c2c2`.
- Settings: `T1_MPRAGE`, four registration levels, 50 registration
  iterations, five post-refresh registrations, registration maximum 5,
  segmentation maximum 50, convergence threshold `1e-4`, and maximum
  distance `1e-4`.

The Java harness ran under a 12 GB watchdog and peaked at **3.504 GB** sampled
RSS. The Python run ran under a 12 GB watchdog and peaked at **8.025 GB**.
Both exited normally. Subject images, raw arrays, classifications, and logs
remain in ignored `build/` and are not committed.

## Reproduction artifacts

- Python run and log:
  `build/invivo-100260-corrected-second-toads-2026-10-14-retry/`
- Java input, source output, and log:
  `build/invivo-100260-source-reference-2026-10-14/`
- Machine-readable comparison:
  [comparison metrics](in-vivo-100260-second-toads-java-comparison-2026-10-14.json)
- Full-suite evidence for the unchanged implementation:
  [LPS compatibility suite report](full-pytest-2026-10-14-toads-lps-axis-compatibility.md)
