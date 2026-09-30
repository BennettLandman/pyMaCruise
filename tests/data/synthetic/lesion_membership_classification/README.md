# LesionTOADS membership classification phantom

This synthetic asymmetric 3 x 2 x 2 phantom pins the selected VandyCRUISE
`LesionToadSegmentation.updateLesionClassifcationformMemberships()` excerpt
(source SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`).
The Java harness is an independently compiled excerpt, not a full
LesionTOADS run. It covers strict equality at the decision boundary, signed
classification labels, an absent Cerebellum-GM name (Java default index 0),
multiple case-insensitive matches (last match wins), and `hasLesion == false`.

The three stored input images are NIfTI-GZip with identity RAS affine and
one-millimeter voxel axes: `classification.nii.gz` is signed byte, while
`memberships.nii.gz` is float32 with the class axis last and
`lesion_membership.nii.gz` is float32. The reference JSON records the exact
inputs and each harness output. Regenerate it with
`sh scripts/generate_java_lesion_membership_classification_reference.sh`, then
regenerate image fixtures with `python tests/data/synthetic/lesion_membership_classification/generate_fixtures.py`.

This method produces an available membership-derived binary short map. No
`PlugInLesionToads` caller was found in the supplied selected source tree, so
this result does not imply the map participates in the active application
workflow. It does not implement the larger LesionTOADS segmentation.
