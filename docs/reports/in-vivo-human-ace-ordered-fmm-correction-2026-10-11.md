# Human MaACE source-ordered FMM correction — 2026-10-11

## Cause and correction

The pre-fix human trace first diverged at Java `02_PHIinit`: 602,421 of
3,195,483 distance values differed, with MAE 2.1245e-5 and maximum error
0.0282317 mm. Java `DistanceFieldyh` seeds and updates an indexed
`BinaryMinHeap` of float-priority voxels in k/j/i scan order and the
`MaskVolume6` neighbor order. The former Python implementation used a
duplicate-entry `heapq` with array-coordinate tuple tie-breaking and did not
preserve the Java insertion/update order.

`_signed_distance_fmm` now delegates to the already source-ordered Numba unit
FMM used by TGDM. This retains the Java float-priority indexed-heap behavior
and its scan order. It also explicitly returns the signed distance cap when an
input has no zero or sign-crossing seed, matching `DistanceFieldyh`'s
FARAWAY-to-cap finalization. The shared operation is used by ACE, ACE narrow
band reinitialization, ThinACE, and Dura.

On the saved human fixed input, the corrected initial level set matches Java
`01_tmpVolfirst` exactly, and the corrected initial distance matches Java
`02_PHIinit` exactly (0/3,195,483 differences; MAE and max error 0). Replaying
ACE on the same fixed human inputs then matched the saved Java outputs at every
voxel:

| Product | Changed voxels | MAE | Maximum error | Dice |
|---|---:|---:|---:|---:|
| Enhanced GM | 0 / 3,195,483 | 0 | 0 | — |
| Cleaned skeleton | 0 / 3,195,483 | 0 | 0 | 1.0 |
| Thinned skeleton | 0 / 3,195,483 | 0 | 0 | 1.0 |

This is a held-fixed ACE source comparison on real human MRI-derived inputs. It
does not rerun Java, TOADS, tissue preparation, or the complete workflow, and
is not clinical validation.

## Validation and resources

The focused affected ACE, Dura, skeleton, and shared unit-FMM checks passed:
**66 passed, 2 expected RAS-orientation warnings, 22.56 seconds**. The run
used the 16 GB RSS watchdog and peaked at **0.363 GB**. This included
`test_ace_reference.py`, whose local saved-human regression compares both the
initial level-set input and its signed-distance output exactly when the
ignored saved human/Java build artifacts are available; that test skips when
those local artifacts are absent. It also included the committed synthetic
Java ACE stage references, ACE after-Gaussian case, ThinACE references, Dura
references, and unit FMM tests.

The fixed-input Python ACE replay command was:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python \
  build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/stage-localization/rerun_fixed.py
```

It completed in **120.61 seconds**, peaked at **0.809 GB**, and produced the
exact metrics above. Detailed JSON and generated Python outputs are in the
ignored `build/invivo-base-source-v3-human-boundary-2026-10-11/ace-java-held-fixed/stage-localization/`
directory. The focused test command was:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/pytest -q \
  tests/test_ace.py tests/test_ace_reference.py \
  tests/test_ace_after_gaussian_reference.py tests/test_skeleton_reference.py \
  tests/test_toads_dura.py tests/test_toads_dura_outputs.py \
  tests/test_dura_wrapper_41.py tests/test_dura_folded_cortex_41.py \
  tests/test_plus_fine_grid_dura_java_reference.py tests/test_tgdm_fmm.py
```

## Source and change provenance

Selected Java reference source SHA-256:

| Source | SHA-256 |
|---|---|
| `DistanceFieldyh.java` | `a0df8859367a52a6273c488fed74d3fbdf0bbf5776b6e6cfea15261a5b290811` |
| `BinaryMinHeap.java` | `5142fbe986301bfbd752c399b98f1997e96e6ec3b1bdbfcc0d0dca14ec472bb5` |
| `MaAnatomicallyConsistentEnhance.java` | `f3cd1c9f4e90cd58ae5cdf37237d90f5a57c59f1a6d6870c9106f4850174cffb` |

The source comparison is localized in
[`DistanceFieldyh.java`](../../tests/java_reference/upstream/edu/vanderbilt/masi/algorithms/CRUISE/ace/DistanceFieldyh.java#L114),
[`BinaryMinHeap.java`](../../tests/java_reference/upstream/edu/jhu/ece/iacl/jist/structures/data/BinaryMinHeap.java#L54),
and [`tgdm_fmm.py`](../../src/pymacruise/tgdm_fmm.py#L287). The Python changes
are in `src/pymacruise/ace.py`, `src/pymacruise/skeleton.py`,
`src/pymacruise/toads_dura.py`, and `tests/test_ace_reference.py`. The test
adds both a small committed Java-reference regression and a local saved-human
stage regression without adding the human images or ignored build products to
the repository.

No full workflow, full suite, or Java rerun was performed. No README or plan
file was edited and no commit was created. The working-tree base commit was
`ecc9c3e`.
