# Atlas-aligned synthetic plus upstream comparisons — 2026-10-10

## Tissue fusion

On the exact 60×75×52 RAS crop from the P41z terminal composition, preserved
`CombineMAandTOADS.java` and Python tissue fusion match exactly on all 234,000
voxels for fused GM, WM, and brain mask. This is selected-source equality for
the tissue-fusion stage; it does not validate TOADS/LesionTOADS segmentation or
registration.

## Topology correction and MaWM refinement

The first comparison showed 14 differences because its Java test harness used
regularization `0.05`. The MaCRUISE wrapper explicitly selects
`1.0e-4`; after correcting the harness setting, the selected Java and Python
topology-corrected field, MaWM erosion, and MaWM-refined field each match
exactly on all 234,000 voxels. The final comparison is saved in the
machine-readable topology report.

The selected Java TGDM wrapper was also run from the Java-refined membership.
Java/Python TGDM surface gates remain unchanged and pass: Dice 1.0, mean
surface distances below 0.00008 mm, and Hausdorff below 0.101 mm for central,
inner, and outer. The algorithm ports match on this crop; the 0.47 mm ribbon
therefore reflects the synthetic inputs and the registered LesionTOADS result
that reached its configured iteration limit without convergence.

## Interpretation

Together with the held-fixed TGDM comparison, these checks show exact tissue
fusion and topology/refinement plus Java-equivalent surface geometry on the
selected crop. They narrow the failed 5 mm ribbon support/quality result to
the upstream synthetic anatomy or registered LesionTOADS output. The registered
loop stopped at its iteration cap, so its output is not a naturally converged
lesion segmentation.

## Evidence and reproduction

- Tissue comparison: [machine-readable report](plus-aligned-folded-tissue-java-2026-10-10.json).
- Topology comparison: [machine-readable report](plus-aligned-folded-topology-java-2026-10-10.json).
- TGDM from the saved Python seed: [report](plus-aligned-folded-tgdm-java-2026-10-10.md).
- TGDM from the Java-refined seed: [machine-readable sensitivity report](plus-aligned-folded-java-topology-seed-sensitivity-2026-10-10.json).
- Ribbon criteria remain failed as documented in the [terminal ribbon audit](plus-aligned-folded-ribbon-quality-2026-10-10.md).

All comparisons are synthetic. They do not establish whole-plus equivalence or
real-MRI validation.
