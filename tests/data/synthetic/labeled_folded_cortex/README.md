# Labeled folded-cortex synthetic case

This committed, deterministic case supplies explicit MaSeg-like tissue labels
and all TGDM inputs on a single 65 × 65 × 65, 1 mm isotropic RAS grid. It is
synthetic integration evidence. It is not a Java-equivalence result or a real
MRI validation.

## Inputs

| File | Meaning and assumptions |
| --- | --- |
| `t1.nii.gz` | Synthetic scalar T1 contrast: background 0, CSF 35, WM 110, GM 145. |
| `brain_mask.nii.gz` | Binary intracranial containment mask on the T1 grid; 255 inside and 0 outside. It includes the synthetic CSF shell. |
| `tgdm_brain_bound.nii.gz` | Inside-negative/outside-positive unit-sign bound derived from the binary mask for `GenericTGDM`; this records the adapter contract used by the Java harness. |
| `maseg_labels.nii.gz` | Discrete labels: 0 background, 5 sulcal CSF, 44 WM, 100 cortical GM. These values explicitly identify the tissue classes for the fixture. |
| `initial_membership.nii.gz` | 255 within the brain mask and 0 outside; seed for the source-selected 126.9 membership iso-level. |
| `white_matter.nii.gz` | 255 in the analytic WM region, 0 elsewhere. |
| `gray_matter.nii.gz`, `ace_gray_matter.nii.gz` | 255 in the analytic 5 mm cortical ribbon, 0 elsewhere. ACE GM is deliberately the same field for this controlled TGDM case. |
| `ace_skeleton.nii.gz` | Zero everywhere to isolate the nested-surface forces; MaACE and skeleton behavior are covered by separate source comparisons. |
| `ground_truth_inner_levelset.nii.gz`, `ground_truth_outer_levelset.nii.gz` | Analytic signed fields for the label geometry; these are descriptive references, not the evolved-surface target. |

The outer cortical boundary is a smooth, low-order angularly folded surface.
The case has no noise, registration, atlas uncertainty, or nonzero GVF; the
zero GVF is explicitly selected for the controlled TGDM run. See `spec.json`
for the analytic shape, source settings, and SHA-256 inventory.

## Run and acceptance checks

From the repository root:

```sh
./.venv/bin/python scripts/run_labeled_folded_cortex.py \
  --output-dir build/labeled-folded-cortex-source-selected
```

The command runs the complete selected source iteration profile, with WM iso
0.6 as selected by the `run_ToadsCruise` wrapper, and writes RAS inner,
central, and outer level sets plus `ribbon-quality.json`. It checks finite
same-grid inputs, label presence, closed/manifold single-component surfaces,
Euler characteristic 2, voxelwise inner ⊆ central ⊆ outer nesting, and finite
positive symmetric nearest-vertex distances between the inner and outer
surfaces. It also requires at least 97% of WM labels inside the inner surface,
GM labels between the inner and outer surfaces, and CSF labels outside the
outer surface. The approximate symmetric nearest-vertex thickness must have a
median within 1 mm of the analytic 5 mm ribbon, with its 5th/95th percentiles
between 2 and 8 mm. The report records exact disagreement counts and thickness
percentiles so they remain visible rather than hidden by a pass flag.
The vertex-distance thickness is an approximate diagnostic, not exact
point-to-triangle cortical thickness.

The agreed source-comparison gates (per-region Dice ≥ 0.97, symmetric mean
surface distance ≤ 0.01 mm, and Hausdorff distance ≤ 1 mm) are recorded but
are not evaluated by the Python-only command. Selected Java outputs for the
same case are committed under `java_reference/`; regenerate them with
`scripts/run_labeled_folded_cortex_java.sh`, then compare stages with:

```sh
./.venv/bin/python scripts/compare_labeled_folded_cortex_java.py \
  build/labeled-folded-cortex-source-selected \
  tests/data/synthetic/labeled_folded_cortex/java_reference \
  --report build/labeled-folded-cortex-source-selected/java-comparison.json
```

Synthetic quality checks, source equivalence, and real-data validation are
separate evidence categories. The current selected Java comparison passes the
Dice and Hausdorff limits, but its final-surface mean distances are above
0.01 mm; the report records the numerical mismatch without relaxing that gate.
