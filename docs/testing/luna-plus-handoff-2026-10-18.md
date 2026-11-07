# Luna assignment: finish the independently actionable PLUS work

> Historical dated record. Current testing readiness, native/Python atlas preparation and run instructions supersede stale open-work statements below; see [the October 27 readiness report](../reports/plus-testing-readiness-2026-10-27.md).

> **Owner correction, 2026-10-18:** PLUS targets `ibsr11objV3`; the earlier `ibsr8obj` requirement below is superseded. Native assets are present, but descriptor adaptation, Python NIfTI conversion, and matched validation remain open. The existing runtime commands still reflect the historical configuration and need updating after compatibility work. See [selection findings](../reports/plus-atlas-selection-2026-10-18.md).

Prepared October 18, 2026. Treat this as a continuation of existing work.

## Objective and authorization

> **Owner correction, 2026-10-18:** the PLUS AtlasEM target is now
> `ibsr11objV3`, replacing the earlier `ibsr8obj` requirement throughout this
> handoff. Its native source files are present and packaged; checked
> descriptor adaptation, coordinate conversion, and selected-atlas validation
> are still required. The historical wrapper source still names `ibsr8obj`,
> so reconcile that runtime selection before claiming a full comparison. See
> [atlas selection findings](../reports/plus-atlas-selection-2026-10-18.md).

Complete the remaining PyMaCRUISE PLUS implementation, source comparisons,
validation tooling, and documentation that can proceed without the missing
selected AtlasEM atlas package or an unavailable complete original runtime.
Do the work, rather than stopping at an audit or proposed plan. Continue past
individual external blockers and finish the independent tasks. Do not claim
that PLUS is complete until its stated acceptance requirements are met.

Work in `/Users/landmaba/dev/mipav/pymacruise`. Read applicable `AGENTS.md`
instructions and the current checkout before editing. Preserve existing
uncommitted work. Do not reset or discard it. Do not spawn other agents,
push, publish, message other people, or modify Codex memory. Leave the final
reviewed changes in the working tree; report what remains uncommitted.

The previous Luna agent was stopped so the user could start this assignment
manually. Its earlier attempt hit a usage limit. Some tool operations may
have been interrupted; establish actual test/process state before launching
another expensive run. Do not stop unrelated processes.

## Current state: inspect and verify, do not assume completion

The last committed HEAD inspected was `dad77581` (`Record BASE TGDM
mask-level-set validation`). The working tree contains PLUS changes from
the earlier attempt. There is no completed full-suite receipt for these
uncommitted changes.

The committed implementation already supports an optional AtlasEM → Dura
stage before LesionTOADS. A generated-atlas synthetic CLI run completed and
resumed all ten then-existing stages. This is integration evidence, not
selected-atlas or whole-workflow equivalence. Do not repeat the older claim
that AtlasEM and pre-LesionTOADS Dura are entirely absent.

The uncommitted work adds final MaSeg refinement, making the baseline route
ten stages and the AtlasEM-enabled route eleven stages. These counts must
be verified against the final implementation. Historical ten-stage test
receipts describe their historical code; do not retroactively change their
results to eleven stages.

The latest committed full-suite receipt predating these edits was on
`8438185d49fcb8ebd49626dc249d763b89a02ffc`: **1,683 passed, 4 skipped,
0 failed**, 30 expected RAS reorientation warnings, 1,000.52 seconds, peak
sampled RSS 9.435 GB under a 16 GB watchdog. It is not validation of the
current working tree.

Earlier Luna reported four focused tests passing before its subsequent
edits, and one file-workflow test passing separately. Do not report those as
validation of the final diff. Recheck the actual final tests and logs.

### Existing edits to review and finish

- `src/pymacruise/plus_segmentation_refinement.py` — new restartable final
  `SEG_refine.nii.gz` checkpoint, input/output hashes, implementation hashes,
  rule metrics, failure manifest, and output repair on resume.
- `src/pymacruise/plus_full_workflow.py` — appends `maseg_refinement`, adds
  `SEG_refine` to final outputs and returns `refined_segmentation_path`.
- `src/pymacruise/plus_file_workflow.py` — removes an unsupported TGDM
  containment-mask argument. Check cache invalidation carefully.
- `src/pymacruise/plus_topology_composition.py` — documents optional explicit
  bounds separately from the selected full-wrapper route.
- `src/pymacruise/cli.py` — newly prints the refined segmentation path;
  verify CLI and config-driven behavior.
- `scripts/prepare_plus_tgdm_java_inputs.py` and
  `scripts/compare_plus_tgdm_saved_java.sh` — revised selected-wrapper mask
  semantics; inspect the complete diff and script compatibility.
- `tests/test_plus_segmentation_refinement.py` — now uses the committed
  MATLAB reference fixture, plus resume/corruption and failed-grid cases.
- `tests/test_plus_file_workflow.py` and `tests/test_plus_full_workflow.py`
  — modified handoff, stage-order, final-output, and integration tests.
- `docs/cli-help.md`, `docs/lesion-toads-port-plan.md`, `plan.md`, and
  `docs/reports/in-vivo-plus-tgdm-mask-aware-java-parity-2026-10-17.md`.
- New report:
  `docs/reports/plus-final-refinement-and-tgdm-mask-routing-2026-10-17.md`.

Use `git status --short` and `git diff` to discover any additional changes.

## Source contract and evidence rules

Read `README.md`, `plan.md`, `docs/plus-baseline.md`,
`docs/validation-tolerances.md`, and the relevant sections of
`docs/lesion-toads-port-plan.md`. The selected full profile is
`run_maCRUISEplus_v1.m`, with shared `pipeline/maCRUISEpluspipeline.m`, from
the selected MaCRUISE SVN r9360 source baseline. Preserve its scientific
behavior, including known defects, IEEE behavior, and stopping rules.

The supplied brain mask intentionally replaces BET/SPECTRE. Visualization
and FreeSurfer exports are excluded. Do not add deferred alternative
pipelines or improve anatomy beyond the selected source behavior.

Keep these evidence categories separate: Python integration, selected-source
stage agreement, human execution, and complete independent original-versus-
Python workflow agreement. One category does not establish another.

Required final segmentation agreement is Dice ≥0.97 for every required
label. Each final surface requires symmetric mean nearest-triangle distance
≤0.01 mm and sampled symmetric Hausdorff ≤1 mm in the same physical RAS
frame. Deterministic stages require exact equality where specified. No
universal scalar/thickness tolerance is established: report errors and
unresolved limits; do not invent acceptance thresholds after seeing results.
Anatomical diagnostics do not replace source agreement.

## Work to complete

### 1. Finish and validate final MaSeg refinement

The original full PLUS runner calls `run_SegmentationRefine` after surface
production. The prior Python full runner omitted this final integration.
Review the new checkpoint against the owning MATLAB source and the existing
Python `segmentation_refine.py`; do not duplicate refinement if another
route already performs it.

Verify the selected coefficients, label rules, central filling, iteration
and neighborhood behavior. Preserve and document any existing source
boundary-indexing or early-stop deviations rather than assuming
`source_faithful=True` proves complete fidelity.

Use `tests/data/synthetic/matlab_refinement_reference/` to check actual output
labels against `matlab-refined.nii.gz`, not merely image shape or manifest
presence. Verify final full-grid geometry, exposed parent output and result
path, CLI/config behavior, successful resume, corrupted-output repair, input
or implementation changes, and a useful failed child manifest. Ensure
implementation changes cannot reuse stale outputs. The parent hashes all
Python source files; inspect standalone checkpoint behavior too.

### 2. Close the selected PLUS TGDM mask discrepancy

Important finding: the selected MATLAB surface route is
`run_CruisePlusSurf_pipeline.m` → `run_TGDMinvcrop.m` →
`run_TGDMinvcrop_jist.m` → `edu.vanderbilt.masi.plugins.CRUISE.tgdm.MaTGDMcrop`.
The JIST command supplies initial membership, WM, GM, ACE, skeleton, GVF,
and crop inputs. It supplies **no optional Brain Mask Level Set**.
`MaTGDMcrop.java` initializes that argument to Java `null` when absent and
passes null into its inner, central, and outer solves.

The prior Python full PLUS route injected a binary-mask-derived unit signed
bound. The uncommitted fix omits that TGDM bound for this selected route.
Retain the supplied mask in source-supported preprocessing and other uses.
Retain explicitly requested bounds in lower-level APIs. Do not copy BASE's
Dura level-set routing into PLUS without evidence from its actual wrapper.

Inspect the source chain directly. `MedicAlgorithmNestedTGDM.java` is a
different entry point; its optional-mask behavior alone does not establish
the selected PLUS contract. Newly written MATLAB reference drivers also do
not supersede the original source.

The earlier mask-aware human TGDM receipt remains a valid *explicit-bound,
held-fixed* comparison. It does not validate the selected unbounded wrapper.
Preserve that distinction in the README, plan, reports, and replay scripts.
Run a fresh selected-source Java/Python comparison with the proper unbounded
handoff using existing human or synthetic inputs where feasible. Selected
Java harness execution is available locally even though the complete
MATLAB/MIPAV/JIST application runtime is unavailable.

Check that `plus_file_workflow.py` records the changed TGDM policy in its
settings/dependency signature so direct `--resume` cannot silently reuse
old mask-bounded products. The full runner's implementation hash does not
automatically protect every independently invoked child workflow. Test
numerical handoff/effect and cache behavior meaningfully.

### 3. Audit and finish remaining source-composition work

Build a finite worklist from open PLUS-relevant plan rows and the actual
source. Identify missing code separately from missing validation. Finish
independently testable supported branches in AtlasEM, LesionTOADS, native
transfer/crop, refinement, and final surfaces/thickness.

AtlasEM already has many source-excerpt tests and direct Demon rigid
optimizer comparisons, including nonzero updates, caps 1/2/5, and selected
registration-class masks. Do not duplicate those accomplishments. Prioritize
uncovered composed initialization → registration → recrop → EM/statistics/
prior updates → stopping → output → Dura contracts. Use deterministic
generated fixtures with independent original Java execution where possible;
Python-versus-Python recomposition is only integration evidence.

Investigate discrepancies and fix supported behavior using owning source
evidence. Retain source NaNs and capped/nonconvergent cases when appropriate.
A generated atlas is a test fixture, never a replacement production atlas.
If a task cannot proceed, document its specific dependency and continue the
other tasks instead of labeling all remaining work blocked.

### 4. Check FLAIR registration independently of the missing atlas

The original selected wrapper uses NiftyReg Aladin rigid-only FLAIR-to-T1
registration; Python uses ANTsPy. Determine whether actual Aladin execution
is available. If available, compare a meaningful deterministic case and
record transforms and resampled images. Matching transform direction alone
does not establish optimizer equivalence.

If Aladin is unavailable, finish source-supported transform direction,
reference grid, interpolation, geometry, and provenance checks that can run
locally, and record the exact remaining external dependency. Do not silently
change the selected registration backend or claim equivalence.

### 5. Build the PLUS end-to-end comparison handoff

Use `scripts/compare_external_base_run.py`, its tests,
`docs/testing/full-base-walkthrough.md`,
`docs/testing/end-to-end-external.md`, and `source-reference-run/README.md`
as starting points. Implement reusable PLUS comparison tooling and synthetic
verification without disturbing BASE behavior.

Verify input identity including FLAIR, parent/child hashes, terminal status,
selected AtlasEM stage presence, active output attempts/revisions, grids,
frames, units, required labels, all final surfaces, level sets, and thickness
summaries. Reject missing/corrupt/incomplete or ambiguous products. Preserve
valid failed comparison reports; distinguish tool failure from failed
numerical agreement. Do not assume the BASE container provides PLUS.

Document both independent full runs, the declared input/preprocessing
boundary, exact required atlas/runtime inputs, commands, provenance,
comparison outputs, visual review, and return files. Clearly identify
untested launch/output-layout assumptions. Never seed either scientific
workflow with the other's intermediates and call it an independent run.
Do not present a fully executable selected-atlas PLUS walkthrough until the
required assets/runtime actually exist and the commands are verified.

### 6. Reconcile documentation and record remaining work

Update current README/plan/CLI status and add a dated evidence report for the
final diff. Some October 16 audit/protocol text predates AtlasEM integration;
mark superseded findings clearly while preserving historical evidence.
Current README also predates the uncommitted final-refinement and selected
TGDM findings. Avoid completion percentages unsupported by runtime/source
coverage. Do not change historical results to match new stage counts.

For each work item, report implemented and verified, externally blocked, or
still independently unfinished. Do not hide unfinished implementation under
the missing-file blocker.

## Useful source and evidence locations

- Original MATLAB source is available read-only under
  `/Users/landmaba/dev/deleteme/masimatlab/trunk/xnatspiders/matlab/MaCRUISE/`.
  Verify the selected revision/file hashes before relying on a copy.
- Java sources: sibling `toads-cruise-scm-2026-09-26`,
  `jist-scm-2026-09-26`, `mipav-src`, plus
  `tests/java_reference/upstream` and reference harnesses.
- Existing AtlasEM reports:
  `docs/reports/atlasem-direct-demon-trajectory-2026-10-17.md`,
  `docs/reports/atlasem-direct-demon-pyramid-optimizer-2026-10-17.md`,
  `docs/reports/plus-atlasem-dura-integration-2026-10-17.md`, and
  `docs/reports/plus-atlasem-end-to-end-synthetic-2026-10-17.md`.
- Human held-fixed evidence:
  `docs/reports/in-vivo-plus-tgdm-mask-aware-java-parity-2026-10-17.md`
  and its JSON, with the explicit-bound scope correction above.
- Existing human PLUS outputs are under ignored
  `build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/`;
  TGDM replay data are under ignored `build/human-plus-tgdm-*` directories.
  Inspect manifests and provenance before using cached arrays.

## External blockers that must remain explicit

The selected resource is
`Atlas/ibsr8obj/cruise-atlas-8obj-noBG-EM.txt` **and every referenced payload**,
with attribution/redistribution terms and hashes. It is not just one text
file. Neither a generic atlas nor the bundled LesionTOADS atlas substitutes
for it. Selected-atlas execution and human AtlasEM convergence await it.

Complete independent workflow equivalence also requires a known-working
original PLUS runtime, exact source/version/settings, full logs and
intermediates, and matched input provenance. The standalone Java harnesses
do not replace that full original application run. Record any additional
Aladin/runtime and thickness-correspondence/acceptance dependencies.

## Verification and delivery

Run meaningful focused tests after each implementation group. Relevant
starting tests include `tests/test_plus_segmentation_refinement.py`,
`tests/test_plus_file_workflow.py`, `tests/test_plus_full_workflow.py`, the
affected TGDM replay/preparation tests, and relevant AtlasEM tests. Verify
both baseline and AtlasEM-enabled CLI/config routes and resume/corruption
handling. Reuse committed fixtures where suitable.

Once code is stable, run the full suite once under the existing memory
watchdog. Do not launch concurrent memory-heavy runs. From the repository:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/pytest -q --tb=short
```

Retain logs and report exact tested code identity, command, pass/fail/skip
counts, warnings, runtime, and memory when measured. A failed run or absent
receipt must not become a pass claim. Fix implementation failures and rerun
affected checks; do not repeat the full suite after documentation-only edits.

Finish with a concise report of changes, source evidence, tests, remaining
external blockers, and any independently unfinished tasks. Include clickable
paths to reports and important changes. Do not declare PLUS finished based
on synthetic execution or held-fixed stage agreement.
