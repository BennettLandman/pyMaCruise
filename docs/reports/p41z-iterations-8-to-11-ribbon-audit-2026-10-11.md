# P41z registered continuation and ribbon audit — 2026-10-11

## Result

The P41z synthetic registered-loop checkpoint at iteration 7 was continued for
four passes (iterations 8–11) using the selected VandyCRUISE profile. The run
completed in **12 min 05 s** under an **8 GB RSS watchdog**, with **1.083 GB**
peak sampled RSS. The completed iteration-11 manifest reports energy ratio
`0.1330575943`, well above the `1e-4` convergence threshold; it is a
nonterminal checkpoint.

The saved checkpoint was composed through Dura transfer, MaSeg separation,
and the selected surface stages, then checked with the explicit ribbon audit.
It still fails:

- GM between inner and outer surfaces: **0.1123%** (minimum 97%).
- WM inside inner surface: **59.5109%** (minimum 97%).
- CSF outside outer surface: **99.9635%**.
- Each mesh is watertight, but has **24 connected components** and Euler
  characteristic **14** rather than one connected genus-zero component.
- Closest-point triangle thickness: median **0.3714 mm** (p05 **0.1536 mm**,
  p95 **0.6918 mm**) versus nominal 5 mm.
- Level-set crossing mean thickness: **0.9094 mm**.

Compared with the iteration-6 checkpoint, GM support changes from 0.0990% to
0.1123%, WM support from 59.4499% to 59.5109%, and median triangle thickness
from 0.3741 mm to 0.3714 mm. These small changes do not represent a quality
recovery; the connected-component count rises from 4 to 24. This continuation
is a Python synthetic integration result, not a Java loop comparison, source
workflow equivalence, or real-MRI validation.

## Reproduction

From the repository root, continue iterations 8–11 from the saved iteration-7
registered checkpoint under the memory watchdog:

```sh
env PYTHONPATH=src PYTHONNOUSERSITE=1 \
  MPLCONFIGDIR=build/p41z-case-normalization/mpl-cache-iter8-11 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
    --limit-gb 8 --interval-seconds 1 --report-interval-seconds 30 -- \
    .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
      --registered-checkpoint-manifest \
        build/p41z-case-normalization/iteration-0007/lesion-registered-loop.json \
      --atlas-descriptor \
        tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
      --output-dir build/p41z-case-normalization/iterations-0008-to-0011 \
      --iterations 99 --thinning-limit 1000 --max-passes 4
```

The completed manifest and stage outputs are under the ignored
`build/p41z-case-normalization/iterations-0008-to-0011/` directory. Surface
composition and quality audit outputs are under the corresponding
`iterations-0008-to-0011-surfaces/` directory. The exact ribbon report is
machine-readable at
[`p41z-iterations-8-to-11-ribbon-quality-2026-10-11.json`](p41z-iterations-8-to-11-ribbon-quality-2026-10-11.json).

## Next work

Trace the loss of cortical GM support from the saved lesion memberships and
registered segmentation through native transfer and MaSeg separation. Preserve
the `1e-4` source convergence threshold and current anatomy checks. A source
match for selected TGDM fields does not explain this nonterminal upstream loop
or the failed ribbon.
