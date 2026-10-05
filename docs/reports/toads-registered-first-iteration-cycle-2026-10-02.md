# Registered TOADS first-iteration cycle (P18cybf)

## Ported behavior

Continue P18cybe's selected first degree-zero exact membership update in
`AlgorithmMultipleToads` source order: compute the float32 energy ratio,
register the 11 atlas shape maps at the current post-refresh rigid transform,
apply convergence and iteration checks, grow on the continuing path, refresh
relations/boundary distances/closest labels, transform shape priors at the
new pose, then enter iteration two through initialization and propagation of
thinning.

The adapter derives state that the selected Java constructor initializes
directly: the per-voxel ordering carried out of the first thinning, its
labels/skeleton/volume and availability arrays, descriptor-ordered output labels, constructor label
map, and initialized boundary-distance state. Maximum iteration and relation
limits remain explicit wrapper parameters. The stage kernels were already
ported; this row composes them and adds no newly countable Java method lines.

## Synthetic input

The committed 13³ RAS test images are in
`tests/data/synthetic/toads_registered_first_iteration_cycle/`. T1 values use
`1200*exp(-r²/24)+50*x+25*y+10*z`; the binary mask is the sphere `r² < 25`
centered at voxel `(6,6,6)`. The test uses the committed 11-class atlas
descriptor and a one-level, one-iteration rigid setup. It chooses iteration
one of two and a zero convergence threshold so the selected branch continues.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_first_iteration_cycle.py
```

Result: **1 passed, 0 failed, 0 skipped in 10.97 seconds**.

The test requires an actual rigid registration result and checks the source
stage suffix through the next thinning pass. It recomputes the energy ratio,
grow result from the first thinning's carried ordering, relation maps, transformed shape
priors, and next thinning directly from their existing kernels, comparing the
arrays and queue trace exactly. Component registration and wrapper-control
methods have separate selected-Java reference tests. This is a synthetic
composition check, not a combined Java wrapper run or full TOADS equivalence.

## Source-limit observation

An initial run on the larger 24×23×23 recropped fixture reached the selected
source relation builder's fixed `MAXLB=500` capacity after growth. The Python
port correctly raises when the selected source table would overflow; the test
therefore uses the smaller committed phantom, which stays within the source
limit. This is a source behavior/input limit to watch during later full-workflow
validation, not a tolerance adjustment.


The original zero-ordering description above was incomplete. The persistent Java ordering field is initialized to zero by the constructor, updated by thinning, and then consumed by growth. This state handoff is corrected and validated against the larger synthetic Vandy run in [the ordering-carry report](toads-ordering-carry-2026-10-05.md).
