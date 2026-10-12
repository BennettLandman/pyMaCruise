# P41z iteration-26 checkpoint and source comparison — 2026-10-11

## Restartable synthetic pass

Continued iteration 25 for one source-profile pass. The run completed under a
4 GB RSS watchdog at **0.761 GB peak RSS**. Its manifest hash is
`2e3aa31ff0f98425a5aaaaaee4853c458cb98a2ac95687935d9ad8c245287d21`; all
**9/9** outputs verify, and all **8/8** NIfTIs are finite and RAS on the
71×89×67 grid. The checkpoint is nonterminal at energy ratio
**0.1323500723** (threshold `0.0001`).

## Selected Java/source comparison

The iteration-25→26 replay completed under a 4 GB watchdog at **0.911 GB peak
RSS**. The selected Java comparisons match Python and the persisted checkpoint:

- Shape sampling: **5,503,849** values, exact.
- Thinning: **267,315** queue events; **21,410,005-byte** serialization,
  exact.
- Post-thinning statistics and membership/energy: exact float32 bits with
  zero NaNs on both sides.
- Registration: transform, rotation, objective, and all 13 bounds, exact.
- Growth: **253,943** events; **17,544,318-byte** serialization, exact.
- Relation refresh: **130** rows; **14** diffusion passes; **848,444-byte**
  serialization, exact.

The comparison is on a fixed, hash-verified checkpoint and is not a complete
Java application run.

## Ribbon-quality audit and loop trajectory

The iteration-26 surface tail completed at **0.635 GB peak RSS**. The audit
fails tissue support, topology, and thickness-tail gates:

- WM/GM/CSF support is **99.419% / 94.469% / 99.266%**; GM is below 97%.
- All meshes are watertight and connected, with Euler characteristic **−14**
  rather than 2.
- Surface-thickness median/p05/p95 is **5.714 / 0.238 / 11.444 mm**. The
  level-set mean is **4.737 mm**, and nesting passes.

The alternating quality profile repeats: iteration 26 returns near iteration
24, while iteration 25 resembles iteration 21. Membership mean absolute
differences are **0.000101** for iteration 24→26 and **0.000195** for 21→25;
the respective segmentations differ at **173** and **240** of 423,373 voxels.
The loop remains nonterminal and its source energy-ratio criterion does not
pass. This is synthetic evidence, not real-human PLUS validation.
