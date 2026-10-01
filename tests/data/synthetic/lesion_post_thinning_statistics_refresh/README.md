# LesionTOADS post-thinning statistics refresh

This comparison composes the P39 first-thinning synthetic state with the
selected wrapper's following `computeCentroids()` → `computeVariances()` →
`computeChannelWeights()` calls. It reuses the folded 17 × 17 × 17 RAS NIfTI,
the initial membership bits, and the final segmentation from
`lesion_initial_thinning_composition`; no additional image fixture is needed.

`LesionPostThinningStatisticsReferenceHarness.java` evaluates the selected
source method excerpts for the exact composed state. The fixture has 13
classes: 12 OBJECT classes and one ordinary LESION class, with no BLACKHOLE
class. The black-hole membership map and fixed initial centroid are zero. The
lesion PowerTable exponents are the selected atlas values for T1_MPRAGE (0.4)
and FLAIR (1.2). The harness consumes big-endian binary state assembled by
`scripts/generate_lesion_post_thinning_statistics_reference.py` from the
committed source NIfTI plus the already pinned P39 first-thinning Java result.
It emits unsigned raw IEEE-754 float32 bit patterns for direct comparison.

Selected source hashes pinned in `java_reference.json`:

- `LesionToadSegmentation.java`: `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`
- `AlgorithmLesionToads.java`: `556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`
- copied `PowerTable.java`: `2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`

Regenerate the reference with:

```sh
.venv/bin/python scripts/generate_lesion_post_thinning_statistics_reference.py
```

The focused test checks initial memberships and first-thinning segmentation
against their separate Java reference, then compares class/auxiliary
centroids, class/auxiliary standard deviations, and channel-weight float bits
against this Java harness. This validates the composed fixed-state continuation
only. It does not cover later membership updates, registration, the complete
LesionTOADS workflow, or real imaging data.
