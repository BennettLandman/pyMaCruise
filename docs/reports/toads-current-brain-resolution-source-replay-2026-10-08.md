# Selected VandyCRUISE TOADS replay on the current brain-resolution run

## Scope

This is a held-fixed source comparison of `AlgorithmMultipleToads` using the
selected `MedicAlgorithmToads08` settings on the current seven-stage synthetic
base workflow. Java receives the exact Python-produced T1 volume used by each
corresponding Python pass, along with the selected source-order atlas. This
isolates the TOADS classification stage; it does not compare upstream N4 or
Dura production.

The Java adapter is
[`ToadsWorkflowSourceHarness.java`](../../tests/java_reference/ToadsWorkflowSourceHarness.java).
It uses the selected wrapper parameters: smoothing `0.2`, maximum iterations
`50`, maximum difference `1e-4`, atlas prior/range `2.0/0.2`, rigid
registration, pyramid parameters `4/50/5/0`, correction enabled, polynomial
degree `3`, kernel `30`, and 18/6 connectivity. The source parameter call is
at `MedicAlgorithmToads08.java:280` in the supplied VandyCRUISE checkout.

Source checksums:

- `AlgorithmMultipleToads.java`: `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`
- `MultipleToadDeformableAtlas.java`: `e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`
- `MedicAlgorithmToads08.java`: `d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115`
- Selected atlas descriptor used by Java: `build/java-reference/source-vandy-atlas.txt`

The input volume is `(172, 216, 164)` at `0.9375 mm` isotropic spacing. NIfTI
values were serialized as float32, little-endian, x-fastest. Java label bytes
were read back in the same x-fastest order using Fortran-order reshape. The
committed comparator is
[`compare_toads_java_classification.py`](../../scripts/compare_toads_java_classification.py).

## Results

| Pass | Held-fixed T1 SHA-256 (`.nii.gz`) | Python result | Java result | Changed voxels | Minimum per-label Dice | Gate |
|---|---|---|---|---:|---:|---|
| First TOADS | `0d08d8419e31a1034423062e5d26f55e0bd34d3ae551dab127bd99d17359f28a` | 11 iterations; classification SHA-256 `5f44f26ab64cb922e59fc3ba9ca22b56969a213446bc5a98c417403e6788b28d` | 11 iterations; raw labels SHA-256 `83357ca2f01e01276d796c07c8bd984a2983e0c9f7df15390336169aecfbbba6` | 851 / 6,092,928 | 0.99522679 (label 6) | Passes Dice ≥ 0.97 |
| Second TOADS | `f3a19b22fb014a1dde19bc6a92646cfb6b2b9418be6f7198aae1be5952166f49` | 48 iterations; classification SHA-256 `25a4b3c5b760073e9bdf730fd5e0d2a7def66a5e5289949dc12870b4342483f8` | 24 iterations; raw labels SHA-256 `a220505ae8e073e1623cd4a1483b408415eccb83aeeec1592cd4cd067439bb7e` | 36,625 / 6,092,928 | 0.95620434 (label 5 / SulcalCSF) | **Fails Dice ≥ 0.97** |

The first pass is within the agreed segmentation tolerance for every label.
The second pass is not: labels 1, 5, 6, 14, 15, 16, 17, 22, 24, and 25 have
Dice scores of 0.99654, 0.95620, 0.98988, 0.98392, 0.99732, 0.99858,
0.99747, 0.99143, 0.99693, and 0.99980, respectively. Label 18 is exact.
The notably different convergence counts (24 versus 48) and the label-5 Dice
failure make current full-resolution second-pass TOADS parity an open blocker
for the base workflow; synthetic integration completion must not be described
as whole-workflow equivalence.

## Iteration-trace diagnosis

To localize the stopping mismatch, a build-only copy of the selected Java
`AlgorithmMultipleToads` source printed its existing float32 energy ratio
after each iteration. The only source edit was this diagnostic print; the
instrumented replay produced the same terminal label hash as each
uninstrumented replay (`83357ca2…aecfbbba6` for pass one;
`a220505a…439bb7e` for pass two). In the second pass, Python and Java ratios
already differ in the first loop: iteration 1 is `0.525829971` / `0.5230326`,
iteration 2 is `0.0409516431` / `0.0383726`, and iteration 10 is `0.00370484893`
/ `0.003672864`. By iteration 23 they are `0.0041031288` / `0.0014557661`; at
iteration 24 Python reports `0.0023892622`, while Java reports `1.4748772e-5`
and stops. The maximum absolute ratio difference over Java's 24 iterations is
`0.002797371`. This shows the discrepancy is in the evolving objective/state,
not merely a different terminal comparison operator. Pinpointing which
membership, field-correction, or registration update first causes the drift
is the next TOADS integration task.

The same print-only trace on the first pass also shows energy-ratio drift even
though its final labels pass the Dice gate: the mean/max absolute difference
over 11 iterations is `0.00331037` / `0.02710345`. At the terminal iteration,
Python reports `8.55285762e-5` and Java `9.7743876e-5`; both stop below
`1e-4`. The first-pass result therefore passes the requested label tolerance,
but does not establish energy-trajectory equivalence.

## Two-iteration held-fixed replay

The shortest supported Python loop run uses two iterations. On the same
second-pass T1 and atlas, the selected Java run (iteration cap 2) and Python
run (iteration cap 2) completed under 16 GB watchdogs with peak RSS of 3.782
GB and 6.240 GB, respectively. Both ratios differ at each step:

| Iteration | Python energy ratio | Java energy ratio |
|---:|---:|---:|
| 1 | `0.525829971` | `0.5230326` |
| 2 | `0.0409516431` | `0.0383726` |

The capped terminal maps differ at 58,039 of 6,092,928 voxels. Per-label Dice
is below 0.97 for SulcalCSF (label 5, `0.915408`) and CerebellarGM (label 14,
`0.958100`). Java raw label SHA-256 is
`8d18b1c54687510d9674c79ed1d6a581c3bd526bee1107eadab8de5bef81dc91`; Python
NIfTI SHA-256 is
`b08d3bae4f1bd67f3a9fe867db4b50690975ba145b8cfa1e7df4a74391396e15`. The
two-iteration cap is for diagnosis only; neither run naturally converged.

The exact Python diagnostic command was:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- ./.venv/bin/pymacruise toads-loop --t1-n4 build/base-workflow-2026-10-08-natural-registration/stages/dura/attempt-01/stripped_t1.nii.gz --brain-mask build/base-workflow-2026-10-08-natural-registration/stages/dura/attempt-01/brain_mask.nii.gz --output build/current-toads-selected-wrapper-2026-10-08/iteration2-python/classification.nii.gz --atlas-descriptor build/base-workflow-2026-10-08-natural-registration/resources/dura-atlas/cruise-atlas-10obj-toads2008.txt --registration-levels 4 --registration-iterations 50 --post-refresh-registration-iterations 50 --registration-max-iterations 50 --maximum-iterations 2
```

This shows the output tolerance is already exceeded by the second full loop
cycle. The first divergent intermediate operation still needs to be isolated.

The selected wrapper passes `firstIter=50` and `topIter=5`; Java sets later
rigid updates to `topIter`. The Python comparison was repeated with the
source-matched `--post-refresh-registration-iterations 5
--registration-max-iterations 5`, instead of `50/50`. On this two-cycle
input, the Python output hash and both energy ratios were unchanged from the
`50/50` diagnostic (`b08d3bae…91396e15`, `0.525829971`, `0.0409516431`). The
registration iteration override is not the cause of the observed mismatch in
this replay. The `toads-loop` CLI default for post-refresh registration was
corrected to `5`; the selected wrapper schedule is now `4/50/5/5` across the
four Python CLI settings. A focused parser regression passes **1 test in
0.58s**. This CLI correction does not resolve the loop-state mismatch.

As a diagnostic only, replaying Java TOADS with the Java Dura output instead
of the held-fixed Python Dura input is not included in the table. That first
attempt used the Dura adapter's raw serialization in the wrong order; after
correcting the conversion, the Java Dura T1 differed from Python at five
voxels and its independent TOADS replay stopped at iteration 30. It is a
different upstream input and does not replace the held-fixed second-pass
comparison above.

## Reproduction and validation category

The Java harness was compiled against the supplied VandyCRUISE/MIPAV/JIST
classes and dependency jars in this workspace. Its executions used the
16 GB RSS watchdog with `-Xmx12g`; peak sampled RSS was 3.943 GB for the first
pass and 3.674 GB for the second. The calls used the committed harness and
current workflow inputs under `build/`:

```sh
CP="$(cat build/java-reference/toads-selected-wrapper-harness/classpath.txt):build/java-reference/toads-selected-wrapper-harness"
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- /usr/bin/java -Xms512m -Xmx12g -cp "$CP" ToadsWorkflowSourceHarness build/current-toads-first-java-2026-10-08/current-n4-t1-xfast-f32le.raw build/java-reference/source-vandy-atlas.txt 172 216 164 0.9375 0.9375 0.9375 50 build/current-toads-selected-wrapper-2026-10-08/first-classification-java.raw
./.venv/bin/python scripts/compare_toads_java_classification.py build/current-toads-selected-wrapper-2026-10-08/first-classification-java.raw build/base-workflow-2026-10-08-natural-registration/stages/toads-first/attempt-01/classification.nii.gz
```

The second-pass invocation substitutes the `current-dura-stripped-t1` input
and `second-classification-java.raw` output paths shown in the report table;
the same comparator is used against the saved `toads-second` NIfTI.

This is **selected-source equivalence evidence**, not a general synthetic
regression suite or real-data validation. No numerical algorithm changed in
this evidence cycle; the TOADS CLI parameter default and a focused parser
regression changed, so no full `pytest` suite was rerun. The most
recent full-suite evidence remains 1,502 passed, 0 failed, 0 skipped, and 11
expected RAS warnings on implementation tree `5a0adcf`; see the
[full-suite report](full-pytest-2026-10-08-bundled-atlas-startup.md).
