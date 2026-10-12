# P41z iteration-20 selected Java handoff and complete Python pass — 2026-10-11

## Result

A complete one-pass Python continuation from the hash-verified iteration-19
checkpoint produced the iteration-20 registered-loop state. The run used the
public `plus-lesion-registered-loop` command, completed in about three minutes,
and peaked at **0.777 GB RSS** under a **4 GB watchdog**. All **9/9** output
hashes verified. All **8/8** NIfTI outputs are finite and RAS-oriented on the
71×89×67 shared grid. The saved stage order ends with `registerShapes`,
convergence checks, `propagateGrowing`, and `computeRelations`.

The corrected selected Java/Python replay now follows the full source ordering:
it refreshes the Python distance maps after thinning, as the workflow does,
then passes those maps into the selected Java membership reference. It
confirms:

- All **5,503,849** registered shape values match (maximum absolute error
  0.0).
- Standard-13 thinning serialization is byte-identical at **21,406,805
  bytes**, covering **97,532** initial queue entries, **267,251** queue
  events, **3,170** topology-candidate events, and zero remaining queue items.
- All seven post-thinning statistics groups match float32 bits.
- Ventricle, brainstem, and GM distance maps refreshed from the thinned
  segmentation match the complete iteration-20 checkpoint values. The Java
  harness comparison is for the subsequent membership update, not the maps
  themselves.
- Every next-membership map and energy field, including objective and energy
  ratio, matches selected Java float32 bits and the complete Python checkpoint.
  Both Java and Python report zero NaNs.

Iteration 20 has objective **3279.7939453125** and energy ratio
**0.132798091**, above the stop threshold `0.0001`; the loop is nonterminal.
The checkpoint manifest SHA-256 is
`8c58f4ac1ea80e0dec35475e57c23cc3f276b5ad36d79000ce0b2f8570c9b23b`.
The iteration-19 parent manifest SHA-256 is
`b060756ae8d04803c38941bf84a464dfeea59cc00ef831f5c63738ae24c8cab9`.

## Scope and limitations

This is synthetic evidence. The full Python iteration executes registration,
growth, and relation refresh and passes restartable manifest checks, but this
report does **not** claim Java equivalence for those three iteration-20 stages
or for a complete live Java loop. The selected Java replay covers shape
sampling, thinning, statistics, distance refresh behavior, and membership.
This is not real-human validation and does not resolve the separate
ribbon-quality failures.

An initial version of the handoff comparator reused iteration-19 distance
maps and reported a mismatching ratio. The comparator was corrected to refresh
the ventricle, brainstem, and GM maps after thinning, as the production loop
does. The corrected ratio and membership arrays match the completed
iteration-20 checkpoint; the earlier value is superseded.

## Reproduction

Complete restartable Python pass:

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iterations-0016-to-0019/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iteration20-complete \
  --iterations 99 --thinning-limit 1000 --max-passes 1
```

Selected-method Java comparison:

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/p41z-case-normalization/iterations-0016-to-0019/lesion-registered-loop.json \
  --output-dir build/p41z-case-normalization/iteration19-next-thinning-java
```

Both runs verify source checkpoint hashes. Generated files and machine-readable
comparison output are retained in ignored `build/` paths.
