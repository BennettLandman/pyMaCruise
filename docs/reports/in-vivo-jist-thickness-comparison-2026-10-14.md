# Human selected-JIST surface-thickness comparison — 2026-10-14

## Scope

This compares the selected `ClosestPointThicknessForInner` Java source with
PyMaCRUISE on the corrected case-100260 inner, central, and outer surfaces.
The Java source is compiled directly with small JIST-container stubs. Its
meshes come from the same human level sets used in the separate [JIST surface
extraction comparison](in-vivo-jist-surface-extraction-corrected-profile-2026-10-14.md).
This is held-fixed stage validation, not whole-workflow equivalence.

Selected source SHA-256 values: `ClosestPointThicknessForInner.java`
`9fc19e8b9caf3defd8278a7dec1dfd417f37b0431f579da4cd17123d662534b6`,
`SurfaceIntersector.java`
`0450c8cbb76770ef4458f9540980748959e0c6ec1c780f3b638173491f92916d`,
`IntersectorTriangle.java`
`8fee2adf76b6cc74122d1a76914615229fb01fbced3d0ec4694a93907dd54759`,
`KdTree.java` `155a5e0b2a09ed1cb00a73ec99eff2288f75f6cfad55d26e659e4c72db7207f3`,
and `KdTriangle.java`
`66ae0741225bd1f57031d6122866ff5940c5975c09a98e0940a1e1810120faf1`.

## Exact held-fixed result

The Java harness serializes the JIST mesh coordinates as float32. Running the
Python thickness API on those exact Java mesh coordinates reproduces every
tested Java output bit-for-bit: all three inner point fields, all three
central distance/thickness fields, and the nearest-inner-vertex index for all
240,378 central vertices. The selected Java source and Python implementation
therefore agree on the same human mesh inputs.

## Production-mesh sensitivity

The production VTK meshes have the same ordered faces and vertex counts, and
their corresponding vertices differ from the Java meshes by at most
`1.34e-5 mm`. This is well inside the surface-position limits. However, that
small coordinate difference changes nearest-hit choices at a small number of
queries on these highly fragmented meshes (67 components, Euler −198):

| Python output vs Java output | Max absolute difference | Mean absolute difference |
| --- | ---: | ---: |
| Inner → outer distance | 0.000345 mm | 0.00000276 mm |
| Outer → inner return distance | 0.804639 mm | 0.00002808 mm |
| Mean inner thickness | 0.402318 mm | 0.00001490 mm |
| Central projected thickness | 2.281953 mm | 0.001388 mm |
| Central → inner distance | 0.000488 mm | 0.00000260 mm |
| Central → nearest inner vertex distance | 0.110468 mm | 0.000143 mm |
| Central nearest-inner-vertex index | 1,956 / 240,378 differ | — |

The aggregate production and exact-Java-mesh Python thickness means are
`2.65278937 mm` and `2.65279329 mm`; sample standard deviations are
`0.94725857 mm` and `0.94725228 mm`. There is no separate user-approved
thickness-field tolerance. The direct algorithm comparison passes on identical
mesh data; the production-versus-Java point-data differences are retained as
a precision-sensitive mesh interaction, not concealed inside the surface
position gates.

The comparison peaked at `1.145 GB` sampled RSS under an 8 GB watchdog. The
machine-readable values, source hashes, and corresponding-vertex distances
are in ignored `build/invivo-base-v3-dura002-no-regrid-n4off-2026-10-14/java-thickness-comparison/comparison.json`.

## Reproduction

First generate JIST meshes from the fixed human level sets with
[`compare_in_vivo_jist_surfaces.py`](../../scripts/compare_in_vivo_jist_surfaces.py).
Then run:

```bash
PYTHONPATH=src:scripts PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/test-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 5 --report-interval-seconds 60 -- \
  .venv/bin/python scripts/compare_in_vivo_jist_thickness.py \
  --core-dir build/invivo-base-v3-dura002-no-regrid-n4off-2026-10-14/stages/cruise-core/attempt-01 \
  --mesh-dir build/invivo-base-v3-dura002-no-regrid-n4off-2026-10-14/java-surface-comparison \
  --surfaces-dir build/invivo-base-v3-dura002-no-regrid-n4off-2026-10-14/stages/surface-thickness/attempt-01 \
  --build-dir build/invivo-base-v3-dura002-no-regrid-n4off-2026-10-14/java-thickness-comparison \
  --report-path build/invivo-base-v3-dura002-no-regrid-n4off-2026-10-14/java-thickness-comparison/comparison.json
```
