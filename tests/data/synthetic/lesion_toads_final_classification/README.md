# LesionTOADS final classification phantom

This asymmetric 3×2×2 RAS phantom exercises the selected final sequence in
`AlgorithmLesionToads`: transfer segmentation indices through the template
label table, then apply `LesionToadSegmentation.updateLesionClassifcation`.
The test inputs are ordinary `.nii.gz` volumes: `segmentation`, the previous
byte classification, three class memberships, lesion membership, and
black-hole membership. `java_reference.json` contains output arrays and counts
from `LesionToadsFinalClassificationReferenceHarness.java`, a selected-method
excerpt whose source file SHA-256 is recorded in the reference.

The tissue label cases include lesion-over-tissue, black-hole-over-tissue,
equal memberships, lesion/black-hole ties, non-lesion labels, a tissue win,
and `hasLesion=false`. The source uses strict `>` for both the eligibility
test and the `lesion < blackHole` subtype choice, so ties are observable.
The output mask is short-valued with codes 0 (not lesion), 1 (lesion), and 2
(black hole). This fixture validates only final classification transfer and
subtype assignment; it does not validate the upstream membership evolution.
