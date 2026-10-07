# JIST versus Lewiner on plus workflow level sets — 2026-10-06

## Scope

This comparison isolates the mesh-extraction step on the saved final inner,
central, and outer level sets from both naturally converged synthetic plus
runs (36×45×34 and 64×80×64). Each level set is passed unchanged as float32
data to the selected `IsoSurfaceOnGrid.solveOriginal` Java source and to
scikit-image Lewiner. The Java harness checks the exact array hash, dimensions,
and selected settings (connectivity 18, iso-level 0, `inclusive=false`). Its
voxel-index vertices are transformed by the input NIfTI affine, so both meshes
are compared in RAS millimeters.

Selected JIST source: CVS r1.7,
`IsoSurfaceOnGrid.java`, SHA-256
`a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`.

## Results

| Grid | Surface | JIST V/F | Lewiner V/F | Components / open edges | Symmetric mean distance | Hausdorff distance |
|---|---|---:|---:|---:|---:|---:|
| 36×45×34 | Inner | 998/1,992 | 998/1,992 | 1 / 0 | 1.212×10⁻⁷ mm | 7.170×10⁻⁷ mm |
| 36×45×34 | Central | 998/1,992 | 998/1,992 | 1 / 0 | 1.112×10⁻⁷ mm | 7.673×10⁻⁷ mm |
| 36×45×34 | Outer | 998/1,992 | 998/1,992 | 1 / 0 | 1.125×10⁻⁷ mm | 7.673×10⁻⁷ mm |
| 64×80×64 | Inner | 3,724/7,444 | 3,724/7,444 | 1 / 0 | 1.031×10⁻⁷ mm | 8.279×10⁻⁷ mm |
| 64×80×64 | Central | 3,724/7,444 | 3,724/7,444 | 1 / 0 | 1.037×10⁻⁷ mm | 8.365×10⁻⁷ mm |
| 64×80×64 | Outer | 3,724/7,444 | 3,724/7,444 | 1 / 0 | 1.020×10⁻⁷ mm | 8.432×10⁻⁷ mm |

For these six fixed integration outputs, JIST and Lewiner agree in vertex and
face counts, connected-component count, and edge-closure checks. The maximum
measured mean and Hausdorff distances are below 0.000001 mm. This is well
inside the project surface tolerances, but it is sampled vertex-to-triangle
distance, not a continuous triangle-to-triangle Hausdorff certificate. Equal
counts and near-identical positions do not prove identical face connectivity.

This result separates the two observed effects: the false outer boundary was
caused by the uncrop background bug fixed in `2f863b4`; the remaining
36-versus-64 grid displacement is already present in resampled tissue supports
and is not introduced by the Lewiner surface substitution on these outputs.
This does not validate topology/ACE/GVF/TGDM against the complete Java plus
workflow.

## Reproduction

With the completed synthetic plus manifests in `build/`, run:

```bash
PYTHONPATH=src .venv/bin/python scripts/compare_plus_grid_surface_extraction.py
```

The script compiles the selected source and harness with local JIST stubs,
executes Java with a 512 MB heap cap, validates the input hashes, and writes
per-level-set metrics to the ignored file
`build/java-reference/plus-grid-surface-comparison/comparison.json`.
Validation command used here:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_plus_grid_surface_extraction.py
```

Result: exit 0; six Java/Python mesh comparisons completed; peak sampled RSS
**0.127 GB** under the 4 GB watchdog. The harness stubs progress/image
containers and normal generation; it tests the selected extractor method,
not the complete JIST runtime or VTK writer.
