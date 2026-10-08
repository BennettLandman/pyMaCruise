# Current seven-stage base run: JIST surface extraction comparison

## Scope

This compares selected Java JIST `IsoSurfaceOnGrid.solveOriginal(..., 18,
0.0f, false)` with Python scikit-image Lewiner marching cubes on the exact same
three CRUISE level sets from the full seven-stage brain-resolution synthetic
base run. It addresses the surface-extraction substitution at the current run's
outputs. It does not validate upstream CRUISE equivalence, full workflow
equivalence, the continuous Hausdorff distance, or real MRI.

The held-fixed input fields are finite RAS float32 arrays on a
`(172,216,164)` grid with 0.9375 mm spacing. The comparison computed full-mesh
topology and nearest-vertex-mapped face connectivity, then evaluated nearest
triangle distances on 4,096 deterministic, evenly spaced vertices per
direction for each surface.

## Results

| Surface | JIST vertices/faces/components/χ | Lewiner vertices/faces/components/χ | Sampled symmetric mean (mm) | Sampled maximum (mm) | JIST faces mapped |
|---|---:|---:|---:|---:|---:|
| Inner | 116,166 / 232,328 / 1 / 2 | 116,178 / 232,368 / 2 / -6 | 0.000000074 | 0.000173 | 67.94% |
| Central | 121,074 / 242,144 / 1 / 2 | 121,100 / 242,224 / 2 / -12 | 0.000040749 | 0.333133 | 63.79% |
| Outer | 123,828 / 247,652 / 1 / 2 | 123,872 / 247,828 / 4 / -42 | 0.000000348 | 0.002224 | 67.51% |

The sampled means are below the 0.01 mm target and sampled maxima are below
1 mm. These are sampled estimates, not a continuous Hausdorff certificate.
The Lewiner meshes remain topologically different from JIST: Java produces one
closed genus-zero component for each level set, while Python produces 2, 2,
and 4 components, with different Euler characteristics. Only 63.79–67.94% of
JIST faces map under nearest-vertex correspondence. Thus the geometric sample
passes do not resolve the mesh connectivity deviation.

## Provenance and reproduction

The source fields are the current run outputs in ignored
`build/ace-indexed-heap-2026-10-08/base-seven-stage-run/stages/cruise-core/attempt-01/`.
Their float32 C-order SHA-256 values are:

| Field | SHA-256 |
|---|---|
| Inner | `d038c9331d1d9ebdcb88586689612b293afc7a89118e2ce91e12d88a0a891b4b` |
| Central | `ef6575a19ea51662734d474d48c470c107cfde4591f829e06bfc0d2f2251ef4d` |
| Outer | `8a2b7effc99270da1137d14a4ef89758d86783621d7b9248b552f9c87301d96e` |

Selected `IsoSurfaceOnGrid.java` SHA-256 is
`a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`; the
comparison harness SHA-256 is
`973263fa3e46408c165ad4c707688f7849b030e4be97e6fd8c3829c2cf119be2`.
The command below writes raw inputs, Java meshes, and JSON metrics under
ignored `build/`:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_workflow_jist_surfaces.py \
  --core-dir build/ace-indexed-heap-2026-10-08/base-seven-stage-run/stages/cruise-core/attempt-01 \
  --work-dir build/base-seven-stage-jist-surface-compare-2026-10-08
```

The command exited 0, reported all three input hashes as consumed by Java,
and peaked at **0.564 GB sampled RSS** under the 8 GB watchdog. The full JSON
metrics are reproducible at `build/base-seven-stage-jist-surface-compare-2026-10-08/comparison.json`.
No pytest or full suite was run for this comparison; this is direct synthetic
source-comparison evidence and not a regression-suite result.
