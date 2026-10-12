# P41z iteration-25 checkpoint and source comparison — 2026-10-11

## Restartable synthetic pass

Continued the registered LesionTOADS case from iteration 24 for one bounded
pass. The run completed at **0.737 GB peak RSS** under a 4 GB watchdog. The
iteration-25 manifest contains **9/9 verified outputs** and **8/8 finite RAS
NIfTIs** on the 71×89×67 grid. The loop remains nonterminal at ratio
**0.1327985525** (source threshold `0.0001`).

## Selected Java/source comparison

The iteration-24→25 comparator completed at **0.898 GB peak RSS** under a 4 GB
watchdog. Selected Java methods match the Python continuation and persisted
checkpoint:

- Registered shape sampling: **5,503,849** values, exact.
- Thinning: **267,492** queue events; **21,418,855-byte** serialization,
  exact.
- Post-thinning statistics and membership/energy fields: exact float32 bits;
  Java and Python have zero NaNs.
- Registration: objective, transform, rotation, and all **13** shape bounds,
  exact.
- Growth: **254,370** events; **17,561,113-byte** serialization, exact.
- Relation refresh: **130** rows; **14** diffusion passes; **848,444-byte**
  serialization, exact.

This is selected-method evidence on a hash-verified fixed checkpoint, not a
complete Java application run. The comparison scope now correctly identifies
iteration 24 as its input.

## Ribbon-quality audit

The iteration-25 surface tail completed at **0.650 GB peak RSS** and fails the
agreed checks:

- WM/GM/CSF support is **99.430% / 92.937% / 99.258%**; GM is below 97%.
- All meshes are watertight and connected, but Euler characteristic is **−18**
  rather than 2.
- Surface-thickness median/p05/p95 is **6.506 / 0.236 / 12.332 mm**; p05 and
  p95 fail. Level-set mean is **4.846 mm**, and nesting passes.

The near two-cycle persists. Same-phase iteration pairs 22→24 and 21→25 have
membership mean absolute differences **0.000100** and **0.000195**, and only
170 and 240 changed segmentation voxels, respectively. Adjacent states differ
more; iterations 24 and 25 also return to the two alternating ribbon-quality
profiles. This is synthetic evidence of a nonterminal trajectory, not
real-human PLUS validation.
