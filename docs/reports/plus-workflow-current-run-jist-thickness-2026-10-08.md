# Current synthetic plus run: JIST thickness wrapper comparison

## Scope

This runs the selected VandyCRUISE `ClosestPointThicknessForInner` Java
wrapper and the production Python thickness API on the same Java
`IsoSurfaceOnGrid` meshes from the current nine-stage synthetic plus workflow.
It isolates the wrapper arithmetic and triangle-tree behavior with mesh
extraction held fixed; it does not establish whole-workflow source equivalence
or real-MRI behavior.

## Wrapper comparison

Each Java mesh has 3,724 vertices and 7,444 faces. The maximum and mean
absolute error for all seven physical fields are **zero**. The 3,724 central
nearest-inner-vertex indices also match exactly. This passes the selected
wrapper's 0.001 mm physical-field comparison gate.

| Output fields | Maximum absolute difference | Mean absolute difference |
|---|---:|---:|
| Inner-to-outer distance | 0 mm | 0 mm |
| Outer-to-inner return distance | 0 mm | 0 mm |
| Mean inner thickness | 0 mm | 0 mm |
| Central-to-outer distance | 0 mm | 0 mm |
| Central-to-inner distance | 0 mm | 0 mm |
| Central projection distance | 0 mm | 0 mm |
| Central nearest-inner-vertex index | 0 indices | 0 indices |

The Python API took **8.51 s** and peaked at **0.132 GB** process RSS. Java
ran with a 2 GiB heap cap; the enclosing watchdog sampled 0.123 GB peak RSS.

## Effect of mesh substitution

Separately applying the same Python thickness operator to the JIST and Lewiner
mesh sets, extracted from the same current plus level sets, gives:

| Mesh set | Mean thickness (mm) | Standard deviation (mm) |
|---|---:|---:|
| JIST meshes | 0.062394644 | 0.018234614 |
| Lewiner meshes | 0.062409653 | 0.018223307 |
| Lewiner minus JIST | +0.000015009 | -0.000011307 |

This is a geometry-impact measurement using the same Python operator on both
mesh sets; it is not a Java wrapper comparison. P08 does not define a separate
thickness acceptance tolerance. JIST/Python face connectivity still differs
on this fixture despite identical component counts and near-coincident
surface positions; see the [current-run surface report](plus-workflow-current-run-jist-surfaces-2026-10-08.md).

## Provenance and reproduction

The Java reference inputs are the current synthetic plus run's level sets;
their NIfTI SHA-256 values are:

| Field | SHA-256 |
|---|---|
| Inner | `402f90a288a0e4b60bc42ed065c7ebd6da04eadf53f616ab529d01b4d7f92c4c` |
| Central | `9a15b9758221684bd176306ce990d744bf1ba663e67e41d16fb399d60380b41a` |
| Outer | `f7d13a1f9373976f524a3642ac9eae9509026577f145d8413a20a5f830f22d7a` |

The JSON report records selected Java source hashes, mesh hashes, dimensions,
comparison metrics, and runtime. Reproduce after regenerating the current
JIST meshes with `compare_workflow_jist_surfaces.py`:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 4 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_jist_closest_point_thickness_brain_resolution.py \
  --mesh-dir build/plus-current-run-jist-surface-compare-2026-10-08/meshes \
  --levelset-dir build/plus-scaled-grid-current/natural-run-64/code-revisions/b200db7a08b9bec9/plus-surface/stages/topology_ace_gvf_tgdm/attempt-01 \
  --work-dir build/plus-current-run-jist-thickness-compare-2026-10-08 \
  --report-json build/plus-current-run-jist-thickness-compare-2026-10-08/comparison.json \
  --report-markdown build/plus-current-run-jist-thickness-compare-2026-10-08/comparison.md \
  --title "Selected JIST thickness on current synthetic plus meshes"
```

The wrapper comparator exited 0 with `comparison_pass=true`; no pytest or full
suite was run for this comparison.
