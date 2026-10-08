# Current synthetic plus run: JIST surface extraction comparison

## Scope

This compares the selected Java `IsoSurfaceOnGrid.solveOriginal(..., 18,
0.0f, false)` with production Python Lewiner extraction on the exact
inner/central/outer level sets saved by the current nine-stage synthetic plus
run. The run's nested CRUISE core fields are finite RAS arrays on
`(66,76,66)` with anisotropic spacing `(0.3047619, 0.2430380, 0.4571429)` mm.
All 3,724 mesh vertices in each direction were included in the nearest-triangle
distance calculation.

This isolates surface extraction on held-fixed Python level sets. It does not
establish the upstream Java/MATLAB CRUISE source composite, continuous
Hausdorff distance, or real-MRI behavior.

## Results

| Surface | JIST / Lewiner vertices | Components | Euler characteristic | Symmetric mean nearest-surface distance (mm) | Max over all mesh vertices (mm) | JIST faces mapped |
|---|---:|---:|---:|---:|---:|---:|
| Inner | 3,724 / 3,724 | 1 / 1 | 2 / 2 | 0.0000000106 | 0.0000610 | 72.84% |
| Central | 3,724 / 3,724 | 1 / 1 | 2 / 2 | 0.0000000017 | 0.0000010 | 72.14% |
| Outer | 3,724 / 3,724 | 1 / 1 | 2 / 2 | 0.0000000021 | 0.0000019 | 72.14% |

All vertex-to-triangle means are below 0.01 mm and maxima are below 1 mm.
Each pair has the same component count and Euler characteristic. However,
only 72.14–72.84% of JIST faces map under nearest-vertex correspondence, so the
triangulations are not identical. The current plus fixture does not reproduce
the extra Lewiner components seen on the larger brain-resolution base fields;
that base-grid deviation remains documented separately.

## Provenance and reproduction

The current Python level-set arrays' float32 C-order hashes are:

| Field | SHA-256 |
|---|---|
| Inner | `cbbf23e1dd4a5a0a835eebf558b2144d2737346d53df66e5facff45258cddaee` |
| Central | `117e6415fdff9758ce502e399026656ed1015fd1078c76c9f9ba426ea999cf0c` |
| Outer | `6c964f9da6607533732f1b1b6c527d0c696604c95efe5300ae0bfc0c0d988226` |

Selected `IsoSurfaceOnGrid.java` SHA-256 is
`a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`;
the Java harness SHA-256 is
`973263fa3e46408c165ad4c707688f7849b030e4be97e6fd8c3829c2cf119be2`.
Reproduce with the saved current plus run and:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_workflow_jist_surfaces.py \
  --core-dir build/plus-scaled-grid-current/natural-run-64/code-revisions/b200db7a08b9bec9/plus-surface/stages/topology_ace_gvf_tgdm/attempt-01 \
  --work-dir build/plus-current-run-jist-surface-compare-2026-10-08
```

The comparison exited 0 and peaked at **0.133 GB sampled RSS**. Its full
metrics are in the ignored
`build/plus-current-run-jist-surface-compare-2026-10-08/comparison.json`.
No pytest or full suite was run for this source comparison.
