# Compact VandyCRUISE 12-iteration output comparison — 2026-10-06

## Scope

This synthetic check reruns the compact committed VandyCRUISE fixture using the selected 11-structure atlas, 12-iteration cap, one-level atlas registration, one registration iteration, padding 10, and relation limit 3. It compares the final native-grid classification to the selected Java source run. It does not validate real-MRI behavior, the whole MaCRUISE pipeline, or identical intermediate TOADS trajectories.

## Evidence

- Java source: selected VandyCRUISE `AlgorithmMultipleToads.java`, SHA-256 `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d` (same source hash as the pinned two-iteration reference).
- Input hashes: same T1, brain mask, and atlas descriptor hashes as the committed two-iteration Java fixture.
- Java trace log SHA-256: `450ade72d8434c36da4f3a0f446f601ace800c6a5e3b418daff55e5e2ba27666` (local instrumented-run log; its full log remains a build artifact).
- Java final classification SHA-256: `aff69a659d72b3733edf3d6913c8bd7212d445c977dcc4c2dfe2b2cdb7a2e1c3`, with MIPAV x-fastest voxel byte order.
- Focused validation: `./.venv/bin/python -m pytest -q --tb=line tests/test_toads_vandy_iteration_boundary.py` — **2 passed, 0 failed, 0 skipped in 44.03s**.
- Current Python run: same final classification SHA-256 as Java, exact across all 55,080 output voxels.
- Both runs hit the 12-iteration cap; neither claims convergence equivalence.

## Iteration trace

| Iteration | Java objective | Python objective | Java ratio | Python ratio |
|---:|---:|---:|---:|---:|
| 1 | 50.30009 | 50.3000908 | 1.0909222 | 1.09092224 |
| 2 | 43.742954 | 43.7429543 | 0.13944969 | 0.139449686 |
| 3 | 43.892612 | 43.8926086 | 0.0034154672 | 0.00341538037 |
| 4 | 50.072083 | 50.0720901 | 0.1315275 | 0.131527722 |
| 5 | 45.65079 | 45.6507263 | 0.09237691 | 0.0923784748 |
| 6 | 44.462265 | 44.4228401 | 0.026378557 | 0.0272640754 |
| 7 | 43.67956 | 43.6787224 | 0.017760089 | 0.0168922711 |
| 8 | 43.643696 | 43.6437035 | 0.0008214486 | 0.000802060182 |
| 9 | 43.722015 | 43.7206383 | 0.0017929128 | 0.00176124065 |
| 10 | 43.68666 | 43.6865196 | 0.00080894976 | 0.000780683244 |
| 11 | 43.702415 | 43.7088661 | 0.0003605645 | 0.000511388469 |
| 12 | 44.73426 | 46.1163673 | 0.02333523 | 0.0536041185 |

Objectives and ratios are close through iteration 5. They begin to separate at iteration 6 and diverge further by the cap (Java objective 44.73426, Python 46.11637). Despite this trajectory difference, the final class-label arrays are byte-identical on this fixture. This result does not explain the intermediate drift and does not establish matching convergence behavior under another input or cap.

## Reproduction

```bash
./.venv/bin/python -m pytest -q --tb=line tests/test_toads_vandy_iteration_boundary.py
```
