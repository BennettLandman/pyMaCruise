# Atlas-aligned synthetic plus upstream comparisons — 2026-10-10

## Tissue fusion

On the exact 60×75×52 RAS crop from the P41z terminal composition, preserved
`CombineMAandTOADS.java` and Python tissue fusion match exactly on all 234,000
voxels for fused GM, WM, and brain mask. This is selected-source equality for
the tissue-fusion stage; it does not validate TOADS/LesionTOADS segmentation or
registration.

## Topology correction and MaWM refinement

Selected Java/Python topology outputs differ at 14 of 234,000 voxels (maximum
absolute scalar difference 0.801). No voxel changes sign. MaWM erosion matches
exactly; the 14 topology differences carry through to the refined membership.
This is a small but real scalar-output discrepancy, recorded separately from
surface quality.

To determine whether the topology discrepancy affects the observed ribbon,
the selected Java TGDM wrapper was also run from the Java-refined membership
instead of the Python seed. Java/Python TGDM surface gates remain unchanged and
pass: Dice 1.0, mean surface distances below 0.00008 mm, and Hausdorff below
0.101 mm for central, inner, and outer. The topology output difference does
not explain the 0.47 mm ribbon measured in this case.

## Interpretation

Together with the held-fixed TGDM comparison, these checks show exact tissue
fusion and Java-equivalent surface geometry on the selected crop. They narrow
the failed 5 mm ribbon support/quality result to the upstream synthetic
anatomy, the registered LesionTOADS output, or the way those labels are mapped
into tissue priors. The registered loop stopped at its iteration cap, so its
output is not a naturally converged lesion segmentation. The 14-voxel
topology scalar mismatch remains open for separate diagnosis.

## Evidence and reproduction

- Tissue comparison: [machine-readable report](plus-aligned-folded-tissue-java-2026-10-10.json).
- Topology comparison: [machine-readable report](plus-aligned-folded-topology-java-2026-10-10.json).
- TGDM from the saved Python seed: [report](plus-aligned-folded-tgdm-java-2026-10-10.md).
- TGDM from the Java-refined seed: [machine-readable sensitivity report](plus-aligned-folded-java-topology-seed-sensitivity-2026-10-10.json).
- Ribbon criteria remain failed as documented in the [terminal ribbon audit](plus-aligned-folded-ribbon-quality-2026-10-10.md).

All comparisons are synthetic. They do not establish whole-plus equivalence or
real-MRI validation.
