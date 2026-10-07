# Plus surface scaling comparison — 2026-10-06

## Purpose and evidence class

This is a synthetic resolution-scaling check between the naturally converged
36×45×34 and 64×80×64 raw-input plus workflows. Both fixtures are generated
from `tests/data/synthetic/plus_pre_tgdm_composition`; their generators use
nearest-neighbor support resampling and construct affines that preserve the
source physical field of view. Both output mesh headers declare RAS
coordinates in millimeters. This comparison is not Java/MATLAB source
equivalence and uses no human MRI.

## Metric and results

The surfaces were read with `pymacruise.surfaces.read_legacy_vtk` and compared
with `pymacruise.validation.mesh_surface_distance_metrics`. This is the agreed
symmetric mean nearest-surface distance: each mesh's vertices are measured to
the continuous triangles of the other mesh. The same routine reports symmetric
Hausdorff distance. Values are in RAS millimeters.

| Surface | 36-grid vertices | 64-grid vertices | Mean absolute distance | Hausdorff distance | 36→64 mean | 64→36 mean |
|---|---:|---:|---:|---:|---:|---:|
| Inner | 998 | 3,724 | 0.275928 | 0.847921 | 0.247435 | 0.304422 |
| Central | 998 | 3,724 | 0.226665 | 0.742792 | 0.213140 | 0.240189 |
| Outer | 998 | 3,724 | 0.226497 | 0.742792 | 0.213070 | 0.239925 |

The agreed surface-position criteria are mean distance ≤0.01 mm and
Hausdorff distance ≤1 mm. After the correction below, all three surfaces pass
the Hausdorff criterion but their symmetric mean distances remain above 0.01
mm. The directional means are now similar, and the mesh extents lie in the
same physical cortical region. This is consistent with a remaining
resolution/discretization difference, but is not by itself proof of source
behavior or acceptance.

The independently resampled input supports already differ on the two grids.
Using the same RAS-millimeter surface metric on 0.5-isosurfaces of the binary
input masks gives:

| Input support | Mean distance | Hausdorff distance |
|---|---:|---:|
| User brain mask | 0.171817 mm | 0.721856 mm |
| MaSeg label 100 (GM) | 0.169389 mm | 0.721856 mm |
| MaSeg label 4 (WM) | 0.160994 mm | 0.622357 mm |

Both fixture generators independently apply nearest-neighbor resampling from
the same 25³ support masks to their target grids. These values show that
discretization is present before MaCRUISE processing and is of similar scale
to the remaining output mean distance. They do not explain the entire output
difference and do not substitute for a same-grid Java/MATLAB source comparison.

## Diagnosed output-background defect and correction

The first comparison found 11,916 vertices on the 64 grid and large apparent
surface distances. Inspecting the saved zero-level sets showed that the
interior negative regions had the same physical bounds on both grids, while
the 64-grid level sets were exactly zero outside their uncropped region. The
selected surface adapter changes exact zero to `-0.1`; those zeros therefore
created false surfaces along the outer volume boundary. The 36-grid crop
happened to span its full native grid, so it did not expose this defect.

The supplied `MaTGDMcrop.java` source creates its uncrop parameters with
background `20` and applies them to the inner, central, outer, and signed
volumes. `uncrop_plus_tgdm_outputs` now uses this source value for its
uncropped scalar products, while input cropping retains its independent zero
background. The plus-surface tail dependency signature includes this
background so resume invalidates and rebuilds old zero-background meshes.

Focused regression command:

```text
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 --interval-seconds 0.5 --report-interval-seconds 30 -- ./.venv/bin/pytest -q --tb=short tests/test_plus_tgdm_uncrop.py tests/test_plus_post_tgdm.py tests/test_plus_surfaces.py tests/test_plus_file_workflow.py
```

Result: **6 passed, 0 failed, 0 skipped in 19.57s**, peak sampled RSS
**0.200 GB**. A saved-checkpoint `plus-full-run --resume` then rebuilt only the
surface/thickness tail under a 4 GB watchdog, completed with peak sampled RSS
**0.251 GB**, and left the top-level manifest `complete` and terminal. The
regenerated meshes have 3,724 vertices each and are confined to the expected
RAS cortical bounds.

## Interpretation and next step

The 64×80×64 workflow is a successful synthetic integration and resource
checkpoint, but these measurements do **not** establish spatial stability
across grid resolutions. Input-support resampling already contributes
0.161–0.172 mm mean distance; the processed meshes differ by 0.226–0.276 mm.
This quantifies a resolution-sensitivity risk, not source equivalence. Keep the
thresholds fixed and use same-grid Java/MATLAB comparisons for source
acceptance; do not characterize the cross-grid means as a tolerance-level
result.

## Reproduction

The compared outputs are the `inner_surface.vtk`, `central_surface.vtk`, and
`outer_surface.vtk` files under:

```text
build/plus-natural-convergence-complete/plus-surface/stages/surfaces_thickness/attempt-01/
build/plus-scaled-grid-current/natural-run-64/plus-surface/stages/surfaces_thickness/attempt-02/
```

The build outputs are ignored and are not committed. To reproduce the metric
after running the 36-grid natural-convergence fixture and
`./scripts/run_plus_scaled_grid_64_synthetic.sh`, run:

```bash
PYTHONPATH=src ./.venv/bin/python - <<'PY'
from pathlib import Path
from pymacruise.surfaces import read_legacy_vtk
from pymacruise.validation import mesh_surface_distance_metrics

roots = (
    Path("build/plus-natural-convergence-complete/plus-surface/stages/surfaces_thickness/attempt-01"),
    Path("build/plus-scaled-grid-current/natural-run-64/plus-surface/stages/surfaces_thickness/attempt-01"),
)
for name in ("inner", "central", "outer"):
    reference = read_legacy_vtk(roots[0] / f"{name}_surface.vtk")
    candidate = read_legacy_vtk(roots[1] / f"{name}_surface.vtk")
    print(name, mesh_surface_distance_metrics(
        reference.vertices, reference.faces, candidate.vertices, candidate.faces
    ))
PY
```
