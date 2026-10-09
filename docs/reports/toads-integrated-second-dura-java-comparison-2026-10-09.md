# Integrated second-pass TOADS Dura handoff comparison

**Evidence category:** held-fixed selected-source comparison on synthetic workflow arrays; not whole-workflow equivalence or real-data validation.  
**Implementation/harness commit:** `eb04598` (`Add labeled folded cortex TGDM validation case`).

The selected Vandy `AlgorithmMultipleToads`/`MultipleToadSegmentation`
implementation was run on the exact `stripped_t1.nii.gz` Dura handoff input
saved by the corrected 2026-10-09 base workflow. The Java run was capped at 50
iterations, used the selected `DuraRemoval_input` mode, and stayed under a 16
GB RSS watchdog: **3.722 GB peak sampled RSS**. The supplied workflow's second
TOADS pass naturally stopped at iteration 43; this isolated replay compared
the saved Dura input/output boundary and did not rerun the preceding stages.

The Java harness and comparator are committed at
`tests/java_reference/ToadsWorkflowDuraOutputsHarness.java` and
`scripts/compare_toads_dura_outputs.py`. The Java raw outputs use x-fastest
linear order (reshaped in Fortran order). Comparing against
`stages/toads-second/attempt-01/dura-inputs` gives:

| Product | Result |
| --- | --- |
| Classification | Exact, 0 differing voxels of 6,092,928 |
| CSF | Exact after selected wrapper's unconditional ×255 scale |
| GM | Exact after ×255 |
| WM fill | Exact after ×255 |
| WM mask | Exact after ×255 |

The first-pass replay on its own saved Dura input also matched all five
products exactly; see the [first-pass report](toads-integrated-first-dura-java-comparison-2026-10-09.json).
The second-pass full comparator result and SHA-256 values are committed in
[the second-pass JSON report](toads-integrated-second-dura-java-comparison-2026-10-09.json).
The exact little-endian float32 replay input hashes were `5729e75b…` (first
pass) and `40d4717c…` (second pass); both are full-grid held-fixed replay
inputs.
The selected Vandy source hashes are `b6ee1cd1…` for
`AlgorithmMultipleToads.java` and `63e135f3…` for
`MultipleToadSegmentation.java`; the harness hash is `fa099c3a…`. These source
files were only read; no files outside `pymacruise` were changed.

This closes a held-fixed Dura handoff comparison for the corrected synthetic
base run. It does not resolve the separate full-profile TGDM mean-surface
drift reported in the [folded-cortex validation](labeled-folded-cortex-validation-2026-10-09.md).
