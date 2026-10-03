# Full synthetic suite: restartable base workflow checkpoint

- **Command:** `./.venv/bin/pytest -q`
- **Tested implementation commit:** `5b3d1ca69d684514a34b31bbf6f57eafda9b8762`
- **Result:** 1,119 passed, 0 failed, 0 skipped; 11 expected RAS orientation warnings
- **Elapsed time:** 1,365.66 seconds (22:45.66)

This repository-wide run includes the restartable synthetic `base-run` workflow,
its stage-output checks and resume-after-corruption test, along with the TOADS
approximate-membership thinning startup and one selected queue-event transition
added after the prior full-suite snapshot. It also includes earlier folded-
cortex, source-reference, CLI, and geometry regressions.

The base integration is mechanical synthetic evidence: it produces the
documented native-grid RAS `.nii.gz` outputs, stage records, and hashed
`run.json`, and its restart test verifies that unaffected stages are reused.
Both short TOADS invocations in that compact workflow hit their iteration caps
without reaching the source energy threshold. The manifest correctly records
`converged=false`. Neither the integration nor this test-suite run establishes
end-to-end Java equivalence, an accepted segmentation, or real-MRI validation.
The larger folded-cortex test is a Python regression for a composed chain;
its aggregate output is not a single Java-reference comparison.

The suite passed against exactly the commit above; later code or documentation
changes are not covered by this run. No checklist percentage or raw test count
is used as a workflow-completion measure.
