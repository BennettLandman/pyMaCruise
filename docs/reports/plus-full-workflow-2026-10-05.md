# Raw-input-to-terminal plus LesionTOADS workflow

**Implementation commit:** `7dc8a67b10489bbf2ba2e7d743c106366c208ae4`

P41q adds `pymacruise plus-full-run`, a restartable file workflow from T1, FLAIR, a user-provided brain mask matching the T1 grid, and the atlas descriptor through preprocessing, crop, initial atlas alignment, registered recrop, initial LesionTOADS memberships/thinning, the registered branch, and a terminal or continuing registered-loop checkpoint.

The top-level `plus-full-workflow.json` records the ordered five workflow stages, source input hashes (including descriptor-owned topology), settings signature, stage status, child manifest hashes, and terminal state. Each child manifest retains its own output hashes. Resume verifies the same inputs/settings, delegates checkpoint verification and repair to each stage, and rejects changed inputs/settings. The raw FLAIR may be on a different physical grid and is registered to T1. Image outputs are compressed NIfTI in RAS. The required user mask means neither BET nor SPECTRE is invoked.

## Validation

- Focused command: `PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl ./.venv/bin/pytest -q tests/test_plus_full_workflow.py` — **1 passed, 0 failed, 0 skipped in 18.40s** on implementation commit `7dc8a67b10489bbf2ba2e7d743c106366c208ae4`.
- The compact committed synthetic atlas fixture is `tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt`. The test generates T1, shifted-grid FLAIR, and a same-grid mask from committed synthetic NIfTIs, then invokes the CLI through an iteration-limit terminal result.
- Verified stage order and child manifest/output records, T1-then-FLAIR channel order, terminal atlas-grid dimensions and RAS affine, compressed terminal classification/subtype images, label values, valid resume, corrupted terminal-output repair without changing prior child manifests, settings-mismatch rejection, and mask-grid rejection.
- Memory-watched full suite: **1,193 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 995.58s**. Peak sampled RSS was **14.609 GB** under the 24 GB watchdog; the run was not terminated. Command: `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`.

## Scope

This is restartable synthetic Python workflow integration. It is not a selected-Java composite comparison, natural-convergence result on a human-sized image, or real-MRI validation. ANTsPy rigid registration remains a documented backend substitution for NiftyReg/Aladin. The runner currently ends at the LesionTOADS checkpoint; the remaining Dura-to-CRUISE surface workflow and source/output acceptance are still open. No human MRI was used.

Full-suite evidence is in [the matching report](full-pytest-2026-10-05-plus-full-workflow.md).
