# P41z iteration-22 checkpoint and ribbon audit — 2026-10-11

## Restartable synthetic pass

Continued the registered LesionTOADS synthetic fixture from iteration 21 to
iteration 22 using one pass and a 4 GB RSS watchdog. The run completed in about
2.5 minutes at **0.777 GB peak RSS**. The complete manifest hash is
`2dfd92bef7e2f1d748694895a350431a5ccfcb111dd4f3121c25a599f4412624`;
all **9/9** output hashes verify. All **8/8** NIfTI outputs are finite and
RAS-oriented on the 71×89×67 grid.

The loop remains nonterminal at energy ratio **0.1324751973** (threshold
`0.0001`). Iterations 20–22 alternate between objective values near 3,280 and
3,747. Iterations 20 and 22 are much closer to each other than adjacent
iterations: their mean absolute membership difference is **0.000130**, and
only **171/423,373** segmentation voxels differ; their objective differs by
1.514 (0.046%). The adjacent 20-to-21 membership mean absolute difference is
0.003672. This is evidence of a near two-cycle, not proof of its cause; the
source energy-ratio stopping rule compares adjacent objectives and therefore
does not stop on this trajectory.

## Selected Java method comparison

The iteration-21-to-22 replay completed under a 4 GB watchdog at **0.914 GB
peak RSS**. On hash-verified fixed states, selected Java references matched
Python and the persisted iteration-22 checkpoint:

- Shape sampling: **5,503,849** values, exact.
- Thinning: **267,254** events; **21,406,955-byte** serialization, exact.
- Post-thinning statistics and membership/energy fields: exact float32 bits;
  both sides contain zero NaNs.
- Registration transform, rotation, objective, and all 13 shape bounds: exact.
- Growth: **253,868** events; **17,539,403-byte** serialization, exact.
- Relation refresh: **129** rows; **14** diffusion passes; **848,431-byte**
  serialization, exact.

The comparator's report label previously said “iteration 19” for every input.
Its source now formats that claim from the verified checkpoint iteration. The
correction changes report metadata only; the selected method calculations above
were completed on the input and output manifests named here. This is selected
method evidence, not a complete Java application run.

## Ribbon-quality audit

The source-selected surface tail completed at **0.632 GB peak RSS**. It still
fails the agreed quality checks:

- GM-between-surfaces support is **94.429%**, below 97%; WM support is 99.416%
  and CSF support is 99.268%.
- All meshes are watertight and connected, but Euler characteristic is **−16**
  rather than 2.
- Surface-thickness median/p05/p95 is **5.704/0.238/11.438 mm**; p05 and p95
  fail. Level-set mean thickness is **4.735 mm** and nesting passes.

Against iteration 21, GM support and Euler characteristic improve, while the
surface p05/p95 bounds still fail. No real-human acceptance claim follows from
this synthetic run.
