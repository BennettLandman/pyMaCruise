# Human BASE run with default orthogonal RAS regridding

## Result

The restartable BASE workflow completed all eight stages on the supplied
oblique human T1, brain mask, and MaSeg with default regridding enabled. Both
TOADS passes reached their energy thresholds (first: 22 iterations,
`5.49568831e-05`; second: 31 iterations, `8.82653876e-06`) rather than their
iteration limits. The run passed the 16 GB memory watchdog (sampled peak
**7.632 GB**). Its manifest records `status=complete`, and an independent audit
verified all **68 stage-output hashes**; all **51 emitted NIfTIs** were finite
and RAS-oriented.

This is successful real-human workflow integration and default-regridding
evidence. It is not whole-pipeline source equivalence or clinical validation.
The existing same-scan source comparison still fails the agreed segmentation
and surface-position limits. Regridding did not materially change tissue
support from the preceding RAS-grid run, and the run retains many small surface
islands; those remain upstream anatomy/algorithm questions.

## Reproduction

From `pymacruise/`, use the supplied sample inputs and a writable output
directory:

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 MPLCONFIGDIR=build/invivo-regrid-mpl \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/python -m pymacruise.cli base-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/input/T1.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/input/T1_bet_mask.nii.gz \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/input/orig_target_seg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/invivo-base-oblique-default-regrid-2026-10-11 \
  --resume --regrid-inputs \
  --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --crop-edge-width 5 \
  --ace-diffusion-iterations 25 --gvf-iterations 200 --tgdm-iterations 1 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 --wm-iso-level 0.5 \
  --ace-use-updated-csf-in-ace 2
```

The T1 input hash is `aae1e3168beea2f26369d34eea60115675e0f0a0123e9f64df8a00d04cfd068b`,
the mask hash is `c8609e0dee0aa4de9b92f4c102cfd5678290abe1bb18e57415c19cac7c0409ad`,
and the MaSeg hash is `524e088110c616ab261f05f3a09861f2735bde170b5b8495b43118be34005e81`.
The original T1 grid was oblique RAS with shape `184×256×256`. The run generated
a shared orthogonal RAS grid of `207×309×305` at 1-mm spacing. T1 used linear
interpolation; the user brain mask and integer MaSeg used nearest-neighbor
interpolation. The derived inputs and hashes are retained under the ignored
`build/` output and listed in `run.json`.

## Stage and resource checks

| Stage | Result |
|---|---|
| Base tissue preparation | Complete; output checks passed |
| N4 | Complete; output checks passed |
| TOADS first pass | Complete; converged naturally at iteration 22 |
| Dura removal | Complete; output checks passed |
| TOADS second pass | Complete; converged naturally at iteration 31 |
| CRUISE core (ACE/GVF/TGDM) | Complete; output checks passed |
| MaSeg refinement | Complete; source-faithful mode; output checks passed |
| Surface and thickness | Complete; all eight outputs passed checks |

The manifest contains eight completed stages. A second audit recomputed every
stage output SHA-256 and inspected NIfTI finiteness and axis codes; it found no
hash, finite-value, or RAS-orientation issues. The surface/thickness summary
reports mean surface closest-point thickness **2.634 mm** and mean level-set
thickness **4.197 mm**.

## Ribbon quality and held-fixed source check

The human ribbon audit reports nested level sets (`inner` outside `central`:
6 voxels; `central` outside `outer`: 0), GM support between inner and outer of
**82.22%**, and WM support inside inner of **88.38%**. **142 of 143 mesh
components are under 100 vertices** on each of the inner, central, and outer
surfaces. All meshes are closed by edge incidence and have no nonmanifold
edges, but component fragmentation is a material ribbon-quality defect. The
prior aligned-grid BASE run measured 82.30% GM support, 88.46% WM support, and
134 components, so regridding did not materially alter support and the
small-island count changed by nine. These support/topology measures do not
pass the source-output Dice or surface-position acceptance gates.

For the fixed human level sets from this run, a direct selected-JIST
`IsoSurfaceOnGrid.solveOriginal(..., CONNECT_18_6=18, iso=0, inclusive=false)`
replay produced identical vertex/face counts, topology, and connectivity for
all three production meshes. Symmetric vertex-to-triangle mean distances were
`0.00000012`, `0.00000015`, and `0.00000015 mm`; sampled Hausdorff values were
`0.00034527`, `0.00048828`, and `0.00048828 mm`. This passes the selected surface
extractor limits (`0.01 mm` mean and `1 mm` sampled Hausdorff) on fixed fields.
All 143 components occur in both Java and Python extraction, so the islands
are already present in the generated level sets; this comparison does not
explain or approve them.

Machine-readable run manifest, quality audit, and Java comparison are retained
under the ignored `build/invivo-base-oblique-default-regrid-2026-10-11/` tree.
Reproduce the independent audits with:

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 .venv/bin/python \
  scripts/audit_in_vivo_base_products.py \
  build/invivo-base-oblique-default-regrid-2026-10-11/quality-audit-view \
  --maseg build/invivo-base-oblique-default-regrid-2026-10-11/resources/regridded-inputs/MaSeg.nii.gz \
  --brain-mask build/invivo-base-oblique-default-regrid-2026-10-11/resources/regridded-inputs/BrainMask.nii.gz \
  --report build/invivo-base-oblique-default-regrid-2026-10-11/regridding-quality-audit.json

PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.2 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_in_vivo_jist_surfaces.py \
  --core-dir build/invivo-base-oblique-default-regrid-2026-10-11/stages/cruise-core/attempt-01 \
  --surfaces-dir build/invivo-base-oblique-default-regrid-2026-10-11/stages/surface-thickness/attempt-01 \
  --build-dir build/invivo-base-oblique-default-regrid-2026-10-11/jist-surface-comparison \
  --report-path build/invivo-base-oblique-default-regrid-2026-10-11/jist-surface-comparison.md
```

The implementation and integration changes were covered by the full suite on
tested tree `cabb52aa`: **1,554 passed, 4 skipped, 0 failed**, with 11 expected
RAS warnings in 858.73 seconds and 9.443 GB peak sampled RSS under a 16 GB
watchdog. The suite report is
[`full-pytest-2026-10-11-orthogonal-regridding.md`](full-pytest-2026-10-11-orthogonal-regridding.md).
That suite is synthetic/software validation; the human run and fixed-field
Java extraction comparison above are separate evidence categories.
