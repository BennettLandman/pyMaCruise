# Synthetic MaCRUISEplus ribbon-quality audit — 2026-10-09

## Finding

The completed 64×80×64 folded plus workflow is valid synthetic integration and
restart evidence, but this particular fixture is not a useful ribbon-quality
acceptance case. The audit confirms all three output surfaces are closed,
connected, genus zero, and nested. The explicit tissue-label support and ribbon
thickness criteria fail.

The supplied MaSeg image contains labels **0, 4, and 100**, not the dedicated
case's expected background/CSF/WM/GM labels **0, 5, 44, and 100**. Thus the
explicit WM/GM/CSF label-presence check fails. The surface-region support check
also fails: only 0% of the GM-labeled voxels fall between the inner and outer
surfaces (minimum required support is 97%). The final inner, central, and outer
surfaces each have 3,724 vertices and 7,444 faces.
Production triangle closest-point thickness has median **0.0603 mm** (p05
0.0342 mm, p95 0.0964 mm); level-set subtraction has only **7 crossing
voxels** and mean thickness **0.2865 mm**. Both are far below the synthetic
nominal 5 mm ribbon bounds.

This is a fixture-quality failure, not evidence of a Java/Python algorithm
mismatch. The existing plus report compares selected TGDM outputs and JIST
surfaces on held-fixed synthetic inputs; those comparisons do not establish
that this fixture exercises anatomically meaningful ribbon thickness.

## Reproduction

After reproducing the completed run from
[the plus workflow report](plus-workflow-retained-center-2026-10-09.md), run:

```sh
./.venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-retained-center-2026-10-09/plus-full-workflow.json \
  --report build/plus-retained-center-2026-10-09/ribbon-quality-audit.json
```

The audit report is [plus-ribbon-quality-audit-2026-10-09.json](plus-ribbon-quality-audit-2026-10-09.json).
It verifies the input/output hashes, finite RAS grids, closed connected
genus-zero meshes, surface nesting, WM/GM/CSF support within the signed-distance
regions, production triangle thickness, and level-set subtraction thickness.
The current run passes topology, nesting, and finite-grid checks; it fails label
presence, label support, and both thickness checks. The runner records a failed
quality result without changing or suppressing the workflow outputs.

## Next integration milestone

The replacement atlas-aligned fixture is committed at
`tests/data/synthetic/plus_vandy_atlas_aligned/` with labels 0/5/44/100. Its
registered LesionTOADS loop is being continued through the selected 99-iteration
cap. Once terminal, run the persisted Dura/native-transfer/MaSeg/plus-surface
stages and apply these same production ribbon checks. Keep fixture-quality
acceptance separate from selected Java source comparisons, whole-workflow
equivalence, and real-MRI validation.
