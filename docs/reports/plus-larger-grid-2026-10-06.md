# Larger synthetic MaCRUISEplus workflow — 2026-10-06

## Result

The public `plus-full-run` command completed its nine-stage path from raw
synthetic T1/FLAIR/user-mask/MaSeg inputs through terminal LesionTOADS, Dura
outputs and native-grid transfer, cortical separation, and CRUISE surfaces
and thickness on the committed **36×45×34** asymmetric folded-cortex RAS
fixture. The fixture contains 3.5 times as many voxels as the prior 25³
end-to-end plus case and preserves its physical field of view.

The test verifies the four fixture SHA-256 hashes and RAS grids; every
top-level and child stage manifest/output hash; finite RAS `.nii.gz` stage
products; the cropped registered LesionTOADS grid; native T1-grid Dura
products and final level sets/thickness image; the iteration-cap terminal
result; and a successful full-workflow resume. The internal registered crop
is `(32,39,30)`, as expected from topology-support cropping; native Dura and
final CRUISE image outputs are `(36,45,34)` on the original T1 grid.

This is synthetic Python integration evidence. The LesionTOADS loop is
configured to stop at a two-iteration cap; this does not demonstrate natural
convergence, Java/MATLAB whole-workflow equivalence, the agreed numerical
acceptance against the source, brain-resolution robustness, or real-MRI
behavior. No clinical image was used.

## Focused validation

Command:

```bash
env MPLCONFIGDIR=/private/tmp PYTHONPATH=src \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pytest -q --tb=short \
  tests/test_plus_full_workflow.py::test_plus_full_run_completes_on_larger_asymmetric_folded_ras_grid
```

Result on implementation commit `949f4a16acdf590af2f3ad5d52b23b4f19538111`:
**1 passed, 0 failed, 0 skipped in 218.97 seconds**; peak sampled RSS was
**7.511 GB**. The 16 GB watchdog did not terminate the run.

The committed image fixture can be regenerated with
`./.venv/bin/python tests/data/synthetic/plus_larger_grid/generate_fixture.py`.

## Full repository regression

Command:

```bash
env MPLCONFIGDIR=/private/tmp PYTHONPATH=src \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  ./.venv/bin/pytest -q
```

Full suite on tested implementation commit `949f4a16acdf590af2f3ad5d52b23b4f19538111`:
**1,219 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,458.33
seconds**. Peak sampled RSS was **12.697 GB**; the 24 GB watchdog did not
terminate the run. This is Python regression evidence and includes the larger
plus integration test, not whole-workflow source equivalence.

The test adds no selected VandyCRUISE TOADS source lines. The separate
method-level source representation estimate therefore remains **3,511–4,911
of 13,651 lines (about 26–36%)**; this estimate is not runtime Python
coverage or workflow completion.
