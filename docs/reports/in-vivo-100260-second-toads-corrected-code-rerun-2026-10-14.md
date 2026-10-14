# Case 100260 second TOADS pass rerun (2026-10-14)

## Result

The saved 100260 Dura handoff was processed with the current RAS-default
TOADS implementation at its original settings. The registered pass reached its
configured 50-iteration limit without meeting the `1e-4` energy-ratio stop
criterion. The final ratio was `0.00145347673`; the terminal classification
changed 998,708 voxels in the final iteration. The run completed normally and
reported **8.025 GB peak sampled RSS** under a 12 GB watchdog.

The second execution reproduced the earlier 100260 second-pass result exactly.
The classification and all four Dura maps (CSF, GM, WM fill, WM mask) have
identical arrays to the previous run: zero differing voxels and Dice 1.0 for
all 11 labels. Both old and new outputs are RAS on the shared 176×240×256 grid.
This demonstrates deterministic replay of this saved intermediate and confirms
that the recent startup/source-axis change did not resolve its nonconvergence.

This is a focused rerun of one human-image intermediate, not a fresh BASE
workflow. The selected Java comparison on this same saved intermediate is
reported separately in the [TOADS source comparison](in-vivo-100260-second-toads-java-comparison-2026-10-14.md); the available case 100260
source log still does not contain a completed subject-specific full-run output.

## Reproduction

Run from the PyMaCRUISE checkout. Subject inputs and prior outputs remain
outside version control; generated files are ignored under `build/`.

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 1 --report-interval-seconds 60 -- \
  env PYTHONPATH=src ./.venv/bin/python -m pymacruise.cli toads-loop \
  --t1-n4 build/invivo-100260-base-2026-10-13/stages/dura/attempt-01/stripped_t1.nii.gz \
  --brain-mask build/invivo-100260-base-2026-10-13/stages/dura/attempt-01/brain_mask.nii.gz \
  --output build/invivo-100260-corrected-second-toads-2026-10-14-retry/classification-RAS.nii.gz \
  --dura-output-dir build/invivo-100260-corrected-second-toads-2026-10-14-retry/dura-inputs-RAS \
  --atlas-descriptor build/invivo-100260-base-2026-10-13/resources/dura-atlas/cruise-atlas-10obj-toads2008.txt \
  --modality T1_MPRAGE --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --maximum-iterations 50 --convergence-threshold 1e-4 --maximum-distance 1e-4
```

The full progress log and generated NIfTIs are under
`build/invivo-100260-corrected-second-toads-2026-10-14-retry/` and are not
committed. Implementation code was unchanged from tested commit
`0981b8d153b700615450ffa64844b2404d584137`; the full-suite result is recorded
in the [LPS compatibility suite report](full-pytest-2026-10-14-toads-lps-axis-compatibility.md).
