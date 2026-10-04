# Compact source-atlas TOADS loop through its 50-iteration ceiling

## Result

The selected registered Python TOADS loop ran on the committed
`(36,45,34)` source-atlas synthetic fixture using the bundled Vandy descriptor.
Although configured for up to 50 iterations, it met its `1e-4` energy-ratio
threshold and stopped at **iteration 49** with ratio **`1.02321183e-5`**.

This is a successful compact synthetic Python convergence result. It is not a
Java equivalence result: the Java/Python reference comparison currently covers
12 iterations and shows divergence after iteration 2. The run does not validate
brain-resolution convergence, whole-workflow Dice/surface tolerances, or
real-MRI behavior.

Reproduce from the repository root:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --fixture tests/data/synthetic/toads_downsampled_vandy_atlas \
  --atlas data/atlases/cruise-atlas-10obj-toads2008.txt \
  --maximum-iterations 50 \
  --output build/java-reference/python-source-vandy-50iter.nii.gz
```

Registered loop time was **132.165 seconds**. Watchdog peak sampled RSS was
**1.109 GB**, with no termination. The output is finite uint8 RAS NIfTI with
labels `[1,5,6,14,15,16,17,18,22,24,25]` and SHA-256
`ba97a91f4f1be6483644d1a95eecd780c3d0c5eaf3231e9e9996aa584bb5ed14`.
The first 20 energy ratios exactly match the separately repeated 20-iteration
low-memory trace; that run also produced the same output hash in two repeats.
Outputs and logs are ignored build artifacts.

## Interpretation

This result demonstrates that the Python loop can reach its selected stopping
criterion on this deterministic compact synthetic fixture when allowed to run
near the wrapper's iteration ceiling. The earlier brain-resolution runs stop
being comparable by iteration 3, so this does not prove the same behavior on
larger data. Next, locate the post-iteration-2 Java/Python state divergence,
then repeat the stopping/convergence check on the brain-resolution synthetic
base route before claiming P18/P34 acceptance.
