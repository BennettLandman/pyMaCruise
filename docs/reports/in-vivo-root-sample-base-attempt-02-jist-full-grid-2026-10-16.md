# Full-grid human BASE JIST surface and thickness comparison — 2026-10-16

## Result

The production JIST wrapper calls `IsoSurfaceOnGrid.solveOriginal` on each
complete level-set image. Replaying the selected Java method on the full
184×256×256 attempt-02 fields, while retaining Java's float32 point
quantization through the RAS transform, matches the production Python mesh
coordinates and ordered faces exactly for inner, central, and outer surfaces.
All three have zero symmetric surface distance, identical topology (67
components, Euler characteristic −198), and pass the configured 0.01 mm mean
and 1 mm sampled-Hausdorff limits.

The selected `ClosestPointThicknessForInner` Java methods also match the
production Python API bit-for-bit on those full-grid meshes. All seven fields
are exact, including nearest-inner-vertex indices for all 240,378 central
vertices. Thus this held-fixed mesh-stage comparison has no thickness
residual; it does not establish upstream or whole-workflow equivalence.

## Corrected interpretation of the earlier comparison

The earlier harness cropped each level set around its negative voxels. That
kept the same topology but changed float32 interpolation coordinates by as
much as `1.53e-5` voxel. Its comparator also promoted Java's decimal
float32 vertex coordinates to float64 before applying the affine. The
combination produced sparse nearest-hit changes and apparent thickness
outliers up to 2.282 mm. It was a validation-harness coordinate mismatch, not
a production VTK precision defect. The production wrapper operates on full
images; full-grid replay removes the mismatch. The old cropped report remains
as historical evidence and is superseded for source-precision conclusions by
this report.

## Metrics and resources

| Surface | Vertices | Faces | Components | Euler | Mean / sampled Hausdorff |
|---|---:|---:|---:|---:|---:|
| Inner | 240,978 | 482,352 | 67 | −198 | 0 / 0 mm |
| Central | 240,378 | 481,152 | 67 | −198 | 0 / 0 mm |
| Outer | 307,792 | 615,980 | 67 | −198 | 0 / 0 mm |

The full-grid surface replay peaked at **0.984 GB RSS** under an 8 GB watchdog;
the thickness comparison peaked at **1.130 GB** under the same limit. The
mesh comparison retains identical fragmented topology; this is source
agreement, not an anatomical quality claim. Machine-readable thickness
metrics and input/source hashes are in
[`in-vivo-root-sample-base-attempt-02-jist-full-grid-thickness-2026-10-16.json`](in-vivo-root-sample-base-attempt-02-jist-full-grid-thickness-2026-10-16.json).
Focused geometry/thickness tests passed: **18 passed in 6.55 s** with
`tests/test_surfaces.py tests/test_thickness.py
tests/test_jist_closest_point_thickness_report.py`.

## Reproduction

From the repository root, run the full-image selected-source extraction:

```sh
PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/python scripts/compare_in_vivo_jist_surfaces.py \
  --core-dir build/sample-singularity-base-source-profile-2026-10-15/stages/cruise-core/attempt-02 \
  --surfaces-dir build/sample-singularity-base-source-profile-2026-10-15/stages/surface-thickness/attempt-02 \
  --build-dir build/sample-singularity-base-source-profile-2026-10-15/java-surface-full-grid-2026-10-16 \
  --report-path docs/reports/in-vivo-root-sample-base-attempt-02-jist-full-grid-2026-10-16.md
```

Then run the selected thickness source on the generated full-grid meshes:

```sh
PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_in_vivo_jist_thickness.py \
  --core-dir build/sample-singularity-base-source-profile-2026-10-15/stages/cruise-core/attempt-02 \
  --mesh-dir build/sample-singularity-base-source-profile-2026-10-15/java-surface-full-grid-2026-10-16 \
  --surfaces-dir build/sample-singularity-base-source-profile-2026-10-15/stages/surface-thickness/attempt-02 \
  --build-dir build/sample-singularity-base-source-profile-2026-10-15/java-thickness-full-grid-2026-10-16 \
  --report-path docs/reports/in-vivo-root-sample-base-attempt-02-jist-full-grid-thickness-2026-10-16.json
```

The surface comparator defaults to the production full-image grid. Its
`--crop-for-small-reference` option is retained for bounded diagnostic cases;
cropped output must not be used as the reference for precision-sensitive
thickness comparisons.
