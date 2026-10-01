# MaCRUISEplus lesion-atlas selection reference

This is a small source-derived selection check, not a LesionTOADS workflow.
The Java harness is an independently compiled excerpt of the supplied
`MedicAlgorithmLesionToads` atlas-choice switch and descriptor-type check.

The selected MaCRUISEplus MATLAB wrapper supplies `With Lesion` unconditionally.
The chosen bundled 12-class descriptor puts `OuterVentricular-WM` at
zero-based class-array slot 10, with descriptor label 25 and type `lesion`.
The no-lesion T1+FLAIR descriptor has the same class slot and label but type
`obj`. The plugin's explicit no-lesion branch selects its T1+FLAIR descriptor
when both modalities are supplied; otherwise it selects the T1-only resource
path. The latter historical descriptor is not bundled in pymacruise, so this
fixture checks the route path only for that branch.

Label 202 in `../plus_ace_branch/segmentation_label202_overlook_patch.nii.gz`
is a separate MaCRUISEplus boundary-overlook segmentation label. Its paired
label-200 control belongs to that boundary test; neither value selects or
identifies the LesionTOADS atlas class.

`java_reference.json` records the selected source paths and hashes, relevant
line ranges, harness hash, and expected synthetic route cases. This checks
configuration and atlas metadata only; it does not establish a complete
LesionTOADS or MaCRUISEplus result.
