# Full-suite validation: plus LesionTOADS atlas initialization

- **Implementation commit:** `4cf78d507bc161ac397c92361ebf730743dd2b30`
- **Command:** `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,175 passed, 0 failed, 0 skipped; 11 expected RAS reorientation warnings.
- **Duration:** 884.50 seconds.
- **Peak sampled RSS:** 15.604 GB; the 24 GB watchdog did not terminate the run.

The run includes P41g: the restartable plus LesionTOADS initial-atlas-alignment
checkpoint, its bundled synthetic 11³ T1/FLAIR and topology inputs, and the
preceding preprocessing/crop-to-initial-alignment focused chain. The focused
chain passed 7 tests in 25.02 seconds with 1.295 GB peak sampled RSS. The CLI
help was also checked for its crop image, crop manifest, atlas descriptor,
iteration count, pyramid level, output directory, and resume options.

P41g's memberships and initial transform match direct composition of the
selected Python adapters bit-for-bit on the synthetic fixture. The full pytest
run is regression evidence; it is not a full LesionTOADS execution, composite
Java equivalence, validation against the published method, or real-MRI
validation. The fixture is synthetic and does not contain human imaging data.
