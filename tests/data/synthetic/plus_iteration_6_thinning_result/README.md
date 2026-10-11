# P41z synthetic iteration-6 thinning result

These are the Python outputs from the iteration-6 thinning pass after an
exact byte-for-byte comparison with the selected Java reference harness. They
are preserved as RAS `.nii.gz` inputs for the following statistics
comparison. `result.json` records source/harness/input-fixture hashes, class
volumes and minima, the 98,325-entry initial boundary queue, 267,638 queue
events, and hashes for the final segmentation, skeleton, availability, and
ordering maps. The queue was exhausted.

Regenerate the result through the opt-in capture path after a successful
focused comparison:

```sh
JAVA_TOOL_OPTIONS=-Xmx4g \
PYMACRUISE_CAPTURE_ITERATION6_THINNING_DIR=build/plus_iteration_6_thinning_result \
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/pytest tests/test_plus_iteration6_thinning_java_reference.py -q
```

Copy the captured NIfTIs and `result.json` here only after the Java comparison
passes. This output is synthetic selected-method evidence, not real-MRI
validation.
