# 2 mm folded-cortex TGDM resolution case

This fixture is an exact every-other-voxel decimation of the committed
1 mm `labeled_folded_cortex` synthetic case. Its affine doubles voxel spacing
to 2 mm while preserving the physical field of view. It retains the same
labels, analytic nominal 5 mm ribbon, and source-selected TGDM parameters.
This is a resolution-sensitivity experiment, not a new reference anatomy.

Regenerate the committed NIfTIs and hashes with:

```sh
./.venv/bin/python tests/data/synthetic/tgdm_resolution_sensitivity_2mm/generate_fixture.py
```

Run the selected TGDM profile and production ribbon audit with:

```sh
PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 4 --interval-seconds 0.2 --report-interval-seconds 10 -- \
  .venv/bin/python scripts/run_labeled_folded_cortex.py \
  --fixture "$PWD/tests/data/synthetic/tgdm_resolution_sensitivity_2mm" \
  --output-dir "$PWD/build/tgdm-resolution-sensitivity-2mm"
```

The runner exits nonzero because level-set-subtraction thickness falls below
its 5 mm ±1 mm phantom criterion. Other explicit ribbon checks pass. Metrics
and scope are recorded in
[`tgdm-resolution-sensitivity-2026-10-11.md`](../../../../docs/reports/tgdm-resolution-sensitivity-2026-10-11.md).
