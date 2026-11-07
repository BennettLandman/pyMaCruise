# Human PLUS TGDM Java comparison with the brain-mask handoff — 2026-10-17

## Evidence category and scope

This is a **real-human, held-fixed selected-source TGDM comparison** on case 100260's PLUS crop. It is not an independent original MATLAB/MIPAV workflow run, does not validate TOADS/LesionTOADS/AtlasEM, and does not establish full PLUS workflow equivalence.

The source `MedicAlgorithmNestedTGDM` wrapper accepts and passes its optional `Brain Mask Level Set` to `solveInnerSurfaceFromOuter`, `solveCentralSurface`, and `solveOuterSurface`. This report compares an explicitly supplied-bound mode: the Java harness consumes the same signed containment bound as PyMaCRUISE's binary-mask adapter (inside `-1`, outside `+1`) on the same 140×178×134 RAS grid. It remains valid held-fixed evidence for that mode. A later source audit found that the selected PLUS wrapper uses `MaTGDMcrop`, whose optional bound is not supplied by `run_TGDMinvcrop_jist.m`; that selected route therefore passes null. This report does not establish selected PLUS wrapper parity or validate upstream mask generation. See the [source audit and routing correction](plus-final-refinement-and-tgdm-mask-routing-2026-10-17.md).

The Java harness uses the selected source `MaGenericTGDM` and JIST `DistanceField` copies in `tests/java_reference/upstream`. It consumes the hash-recorded precomputed initial level set, tissue fields, ACE-GM, thin skeleton, three-component GVF, and signed mask bound. The Python evolution uses the same arrays and selected wrapper settings. The RAS affine and input hashes are recorded in the [machine-readable comparison](in-vivo-plus-tgdm-mask-aware-java-parity-2026-10-17.json); raw arrays remain in ignored `build/`.

## Result

All compared zero-level surfaces meet the agreed gates: Dice ≥0.97, symmetric mean nearest-surface distance ≤0.01 mm, and Hausdorff distance ≤1 mm. The evolved masks have zero sign-different voxels and Dice 1.0 across initial inner, central, inner, and outer.

| Surface | Dice | Symmetric mean distance (mm) | Hausdorff (mm) | Sign-different voxels |
|---|---:|---:|---:|---:|
| Smoothed initial | 1.000000 | 0.000000365 | 0.006424 | 0 |
| Initial inner | 1.000000 | 0.000000424 | 0.009846 | 0 |
| Central | 1.000000 | 0.000000255 | 0.005421 | 0 |
| Inner | 1.000000 | 0.000000044 | 0.000349 | 0 |
| Outer | 1.000000 | 0.000016482 | 0.226406 | 0 |

The Java replay peaked at **1.356 GB** under an **8 GB RSS watchdog**; Python evolution peaked at **2.493 GB**, and the mesh comparison at **0.571 GB**. The focused regression selection passed **31 tests in 22.76 seconds** before this evidence run. The full suite on implementation commit `50ad2ea1b25ecde845bc8c399813c9da3297b419` passed **1,681 tests, skipped 4, and had 0 failures** in **991.89 seconds**; there were 30 warnings, all from expected RAS reorientation paths. Peak sampled RSS was **7.758 GB** under a **16 GB watchdog**. See the [full-suite evidence](full-pytest-2026-10-17-mask-aware-tgdm.md).

## Reproduction

The preparation, Java replay, and comparison commands used these same case-100260 workflow inputs and this crop:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/plus-surface \
  --context-manifest build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/plus-full-workflow.json \
  --output-dir build/human-plus-tgdm-mask-corrected-2026-10-17 \
  --use-saved-tgdm-outputs --allow-python-stage-differences

javac -d build/human-plus-tgdm-mask-corrected-2026-10-17/classes \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/PlusTGDMReferenceHarness.java

./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 -- \
  /usr/bin/java -Xms256m -Xmx6g \
  -cp build/human-plus-tgdm-mask-corrected-2026-10-17/classes \
  PlusTGDMReferenceHarness \
  build/human-plus-tgdm-mask-corrected-2026-10-17/arrays.bin \
  build/human-plus-tgdm-mask-corrected-2026-10-17/java \
  precomputed_initial_level_set \
  build/human-plus-tgdm-mask-corrected-2026-10-17/brain-mask.bin

PYTHONPATH=src ./.venv/bin/python scripts/compare_plus_tgdm_java.py \
  --workflow-root build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/plus-surface \
  --reference-dir build/human-plus-tgdm-mask-corrected-2026-10-17 \
  --report build/human-plus-tgdm-mask-corrected-2026-10-17/comparison-with-mask-bound.json
```
