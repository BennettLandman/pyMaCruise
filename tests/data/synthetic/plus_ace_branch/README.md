# MaCRUISEplus boundary-to-ACE branch phantom

This committed synthetic RAS fixture provides fixed upstream products for
testing the plus branch from multi-atlas segmentation boundary construction
through tissue-map combination and ACE. It contains a 24³ 1 mm isotropic
concentric tissue phantom and aligned NIfTI inputs in RAS orientation.

`segmentation_label202_overlook_patch.nii.gz` contains a compact patch coded as
label 202, which the selected plus boundary operation treats as an overlooked
label. `segmentation_label200_control.nii.gz` is identical except that patch
retains label 200. Label 202 is not the LesionTOADS atlas class whose descriptor
type is `lesion`; this phantom checks boundary-label filtering only. Both images
are passed through the same boundary operation and TOADS-WM cleanup.
The atlas and TOADS maps are fixed synthetic upstream products; the test does
not run TOADS, LesionTOADS, N4, or any full plus workflow.

For each branch, the test compares raw/cleaned boundaries, updated GM/WM/mask,
and the final ACE outputs with their public stage APIs, while checking that
the input image arrays remain unchanged. The ordinary base ACE call is also
checked against its previously committed Java reference outputs.

Regenerate fixture volumes with `make_phantom.py` from the repository root's
`pymacruise` directory. The ACE gray/white dimensions and analytic geometry
match `../ace/`; this directory separately commits all images needed by the
composed branch and the label-202/label-200 segmentation variants.
