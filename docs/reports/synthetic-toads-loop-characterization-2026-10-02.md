# Synthetic registered TOADS loop characterization

This report records a reproducible synthetic diagnostic for the registered
TOADS loop used by the base-workflow integration milestone. It does not modify
the loop, adjust source thresholds, demonstrate source-equivalent segmentation,
or validate real imaging data.

## Runnable command

From the repository root, run:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --output build/benchmarks/synthetic-toads-10.nii.gz \
  --maximum-iterations 10 \
  --registration-levels 1 \
  --registration-iterations 1 \
  --loop-registration-iterations 1 \
  --padding 2 \
  --relation-limit 1
```

The script reads only the committed 13×13×13 synthetic T1 and mask in
`tests/data/synthetic/toads_registered_first_iteration_cycle/`, plus the
committed 11-class atlas descriptor and images. It executes the existing
atlas-alignment and registered-loop APIs. It reports per-iteration objective,
energy ratio, exact-membership time, registration/branch time, registration
step count, terminal status, and output path/hash. Existing output paths are
left untouched; choose a new `--output` path for each run.

## Observed run

On 2026-10-02, the command above reported:

```text
elapsed.alignment_seconds=0.408
elapsed.loop_startup_seconds=1.320
elapsed.loop_total_seconds=18.342
iteration.1 objective=30.6944065 ratio=1.35140693 membership=1.425s branch=0.318s
iteration.2 objective=20.6838207 ratio=0.389682025 membership=1.421s branch=0.299s
iteration.3 objective=35.2244225 ratio=0.520159483 membership=1.409s branch=0.310s
iteration.4 objective=21.1386795 ratio=0.499821424 membership=1.455s branch=0.312s
iteration.5 objective=37.5722008 ratio=0.559811771 membership=1.412s branch=0.296s
iteration.6 objective=22.2930775 ratio=0.510450244 membership=1.436s branch=0.297s
iteration.7 objective=40.9038429 ratio=0.588976979 membership=1.400s branch=0.293s
iteration.8 objective=23.8644314 ratio=0.526165366 membership=1.444s branch=0.291s
iteration.9 objective=44.2066879 ratio=0.597676575 membership=1.396s branch=0.288s
iteration.10 objective=25.2528152 ratio=0.545753181 membership=1.430s branch=0.089s
terminal iterations=10 ratio=0.545753181 threshold_reached=false iteration_limit_reached=true
output=build/benchmarks/synthetic-toads-10.nii.gz
sha256=779c64ba73059d716fd65fcbbec3193a2e0bd5701cd687e223c530650979d764
shape=(13, 13, 13) dtype=uint8 orientation=('R','A','S') labels=1..11
```

A second direct API run with relation limit 3 and five maximum registration
iterations took 18.682 seconds for 10 outer iterations and ended at energy
ratio 0.5458928, still above threshold. Increasing registration effort on
this fixture therefore did not remove the cap outcome.

## Runtime profile and stop behavior

An instrumented two-iteration run took 5.688 seconds inside the registered
loop. Its `cProfile` output attributed 3.805 seconds cumulatively to three
default-atlas membership passes, 1.042 seconds to inverse-shape-prior energy,
0.738 seconds to shape-relation factors, and 0.673 seconds to PowerTable
lookups. The rigid optimizer took 0.315 seconds cumulatively. This identifies
membership and prior evaluation as the principal runtime targets if later
profiling on larger synthetic volumes confirms the same pattern. The short
run also includes one-time Python imports, so its startup times should not be
scaled linearly to other systems or volume sizes.

The selected VandyCRUISE wrapper computes the ratio after exact memberships
using `abs(2*(Eprev-energy)/(Eprev+energy))`, calls shape registration, then
stops when the ratio is below `maxDistance` or the iteration count reaches
its limit (`VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java`,
lines 766–788). The Python loop uses the same ratio and stop order. In this
phantom the objective alternates between two ranges and the ratio remains
roughly 0.39–0.60; neither registration settings nor the final iteration show
source-threshold convergence. This is synthetic evidence about this fixture,
not proof that the source algorithm converges on every input.

No threshold relaxation, damping, or other algorithm change is justified by
these observations. For the integration milestone, report this fixture's
iteration-limit terminal state explicitly; a separate synthetic case would
need source-grounded construction and validation before it could demonstrate
energy-threshold convergence. The output and trace are not a Java-versus-
Python numerical comparison. No real MRI was used. No focused pytest or full
suite was run for this diagnostic script.

## Bounded synthetic-convergence probes

To check whether the cap result was specific to the original smooth intensity
field or fuzzy shape priors, the trace script also supports intensity and
loop-time shape-prior variants. With the default one-level, one-registration-
iteration, `1e-4` threshold, the reproducible probes were:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py --output build/benchmarks/probe-constant.nii.gz --maximum-iterations 2 --intensity-variant constant-volume --relation-limit 1
./.venv/bin/python scripts/trace_synthetic_toads_loop.py --output build/benchmarks/probe-constant-mask.nii.gz --maximum-iterations 2 --intensity-variant constant-mask --relation-limit 1
./.venv/bin/python scripts/trace_synthetic_toads_loop.py --output build/benchmarks/probe-linear-x.nii.gz --maximum-iterations 2 --intensity-variant linear-x --relation-limit 1
./.venv/bin/python scripts/trace_synthetic_toads_loop.py --output build/benchmarks/probe-linear-xyz.nii.gz --maximum-iterations 2 --intensity-variant linear-xyz --relation-limit 1
./.venv/bin/python scripts/trace_synthetic_toads_loop.py --output build/benchmarks/probe-original-scaled.nii.gz --maximum-iterations 2 --intensity-variant original-scaled --relation-limit 1
./.venv/bin/python scripts/trace_synthetic_toads_loop.py --output build/benchmarks/probe-uniform-shapes.nii.gz --maximum-iterations 5 --shape-variant uniform-one --relation-limit 1
./.venv/bin/python scripts/trace_synthetic_toads_loop.py --output build/benchmarks/probe-onehot-shapes.nii.gz --maximum-iterations 5 --shape-variant topology-one-hot --relation-limit 1
```

Observed terminal ratios were respectively `0.173250198`, `0.240745962`,
`0.238085672`, and `0.164067119`; all runs stopped at the iteration cap, not
at the energy threshold. A constant image inside the brain mask gave the same
two-cycle objectives and ratios as the constant-volume probe because the
preprocessing stage zeros voxels outside the mask. Separate direct-API probes
of linear-x intensity and a 0.001-scaled original intensity ended at ratios
`0.23590134` and `0.38968217` after two cycles, also without convergence.

The shape-prior variants replace shape arrays in the loop state **after** the
ordinary descriptor-based atlas alignment; they characterize sensitivity of
the iteration, not candidate end-to-end atlas assets. The topology-one-hot
variant's five ratios were `0.749594808`, `0.183845341`, `0.080778733`,
`0.223073974`, and `0.164067119`. The uniform-one variant's ratios were
`1.67613351`, `0.216051474`, `0.28869465`, `0.242915615`, and `0.238085672`.
Neither approached `1e-4`.

These bounded probes did not identify a source-faithful, threshold-convergent
synthetic setup. The restartable base runner may be reported as mechanically
complete with its two TOADS invocations explicitly marked capped and
nonconverged; it must not be presented as convergence-validated. A proper
class-aligned synthetic atlas/input design remains a possible future test,
but these experiments provide no evidence that it will converge or that the
current output is source-equivalent. Do not relax the threshold to force the
status flag.
