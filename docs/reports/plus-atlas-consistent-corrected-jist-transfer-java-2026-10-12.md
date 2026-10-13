# Corrected synthetic PLUS JIST transfer comparison — 2026-10-12

## Result

Selected JIST `ReorientVolume.solve(UNCHANGED, NEAREST_NEIGHBOR)` followed by
the wrapper's reference-header assignment is bit-identical to Python's
`jist-unchanged` transfer for CSF, GM, WM-fill, and WM-mask on the corrected
synthetic case. All four channels have zero changed voxels and zero float32
error. The comparison used the saved RAS Dura products and T1 reference.

The source behavior preserves the Dura sample array dimensions, `(166, 206,
160)`, and assigns the larger T1 reference geometry, `(172, 216, 164)`. It does
not resample voxels or expand the array to the reference dimensions. In
contrast, this project's recommended `physical-affine` mode resamples to the
reference grid. Both modes are explicit; this comparison establishes only
that `jist-unchanged` reproduces the selected JIST method. It does not show
which mode was used for the historical in-vivo examples or establish
whole-PLUS equivalence.

The selected Java implementation is run through the existing headless MIPAV /
JIST harness. Java output is converted from MIPAV's x-fastest order to NIfTI
RAS array order before bitwise comparison.

## Reproduction

```sh
PYTHONNOUSERSITE=1 PYTHONPATH=src .venv/bin/python \
  scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 2 \
  --report-interval-seconds 30 -- .venv/bin/python \
  scripts/compare_plus_dura_jist_transfer.py \
  --workflow-root build/plus-atlas-consistent-corrected-full-2026-10-12 \
  --output-dir build/plus-atlas-consistent-corrected-jist-transfer-java-2026-10-12 \
  --report docs/reports/plus-atlas-consistent-corrected-jist-transfer-java-2026-10-12.json
```

Peak sampled RSS was **0.321 GB** under the 8 GB watchdog. The JSON report
includes shapes, Java extents, per-channel equality, and the reference input
hash.
