# Full pytest checkpoint: N4-to-CRUISE synthetic chain

- Command: `./.venv/bin/pytest -q`
- Result: **1,100 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings**
- Duration: **2,218.92 seconds (36:58.92)**
- Tested implementation commit: `a7c142cead48bc13bd1e96854992d903ea3377bf`
- Evaluated: 2026-10-02

This full repository run includes the connected compact synthetic CLI sequence:
actual ANTsPy N4, first registered TOADS, Dura removal, second registered
TOADS, optional same-grid MaSeg cleanup, base ACE, generalized GVF, nested
TGDM, and native-grid inverse crop. The TOADS test iterations and Dura outer
iterations remain intentionally limited on the small phantom. The result
establishes repository-wide test compatibility at this commit; it does not
establish full TOADS convergence, human-data acceptance, historical SPECTRE
space equivalence, or complete Java-wrapper equivalence.

The 11 warnings are the existing expected input-orientation reorientation
warnings for RAS-normalizing commands and their reference calls. No test was
skipped.

This full-suite checkpoint validates the completed, bounded selected-source
slices P18cyau, P18cyav, P18cyaw, P18cyax, P18cyaz–P18cybm, and P18cybo–P18cybx. Each corresponding
Done row in `plan.md` records this exact test result and implementation commit.
The registered-TOADS method slices have their focused Java-reference reports
linked from those rows. The connected workflow evidence is in the [core
composition report](toads-cruise-core-2026-10-02.md), [two-pass
TOADS-to-core report](toads-dura-second-to-cruise-core-2026-10-02.md),
[MaSeg cleanup report](toads-cruise-cerebellum-integration-2026-10-02.md), and
[N4 integration report](n4-toads-dura-cruise-chain-2026-10-02.md).
