# P41z iteration-19 registered-state statistics comparison — 2026-10-11

## Result

The saved iteration-19 registered state was passed through the selected Java
`computeCentroids`, `computeVariances`, and `computeChannelWeights` excerpt and
the Python production statistics seam. All seven float32 result groups match
bit-for-bit: class, lesion, and black-hole centroids; class, lesion, and
black-hole standard deviations; and channel weights. Both implementations
reported **zero NaNs**. The state is 71×89×67 RAS with T1 and FLAIR channels.
The full comparison ran under a 4 GB watchdog and peaked at **0.224 GB RSS**.

The source profile uses `Identity` channel normalization, outlier distance
0.1, optimized CSF power tables on both modalities, tissue power exponents
of one, and lesion exponents 0.4/1.2. The harness source hashes match the
selected VandyCRUISE `LesionToadSegmentation` methods and
`AlgorithmLesionToads` wrapper.

## Scope

This checks statistics on the saved registered segmentation and membership
state. It does not rerun or compare the next iteration's thinning, membership
update, registration, growth, relation refresh, complete registered loop, or
the entire Java workflow. It therefore rules out a statistics mismatch on
this particular saved state, but does not identify the source of the
nonterminal loop or failed ribbon-quality metrics.

## Reproduction

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_p41z_registered_statistics_java.py \
  build/p41z-case-normalization/iterations-0016-to-0019/lesion-registered-loop.json \
  --output-dir build/p41z-case-normalization/iteration19-statistics-java/input \
  --report build/p41z-case-normalization/iteration19-statistics-java/comparison.json
```

The generated registered-state inputs and machine-readable comparison live
in the ignored build directory. The comparator uses the selected-method Java
harness in `tests/java_reference/` and retains all generated artifacts needed
to inspect its inputs.
