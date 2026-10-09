# Synthetic MaCRUISEplus ribbon-quality audit — 2026-10-09

## Finding

The completed 64×80×64 folded plus workflow is valid synthetic integration and
restart evidence, but this particular fixture is not a useful ribbon-quality
acceptance case. The independent audit confirms all three output surfaces are
closed, connected, and genus zero, while the tissue labels and nested surface
separation fail the explicitly selected synthetic ribbon bounds.

The supplied MaSeg image contains labels **0, 4, and 100**, not the dedicated
case's expected background/CSF/WM/GM labels **0, 5, 44, and 100**. The final
inner, central, and outer surfaces each have 3,724 vertices and 7,444 faces.
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
It checks finite RAS grids, closed connected genus-zero meshes, presence of
explicit WM/GM/CSF labels, production triangle thickness, and level-set
subtraction thickness. The current run passes topology and finite-grid checks
and fails the label and both thickness checks. The runner records a failed
quality result without changing or suppressing the workflow outputs.

## Next integration milestone

Create a labeled folded plus input case that carries meaningful WM, GM, and
sulcal-CSF regions through the selected plus workflow. Require the production
JIST surfaces, thickness products, nesting, label support, and the same
synthetic quality bounds used for the dedicated 65³ TGDM case. Keep this
fixture-quality acceptance separate from selected Java source comparisons,
whole-workflow equivalence, and real-MRI validation.
