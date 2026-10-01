# MaCRUISEplus boundary-to-ACE branch phantom

This committed synthetic RAS fixture provides fixed upstream products for
testing the plus branch from multi-atlas segmentation boundary construction
through tissue-map combination and ACE. It contains a 24³ 1 mm isotropic
concentric tissue phantom and aligned NIfTI inputs in RAS orientation.

`segmentation_lesion.nii.gz` contains a compact patch coded as label 202,
which is in the selected plus boundary's overlooked-label list.
`segmentation_no_lesion.nii.gz` is identical except that patch retains label
200. Both are passed through the same boundary operation and TOADS-WM cleanup.
The atlas and TOADS maps are fixed synthetic upstream products; the test does
not run TOADS, LesionTOADS, N4, or any full plus workflow.

For each branch, the test compares raw/cleaned boundaries, updated GM/WM/mask,
and the final ACE outputs with their public stage APIs, while checking that
the input image arrays remain unchanged. The ordinary base ACE call is also
checked against its previously committed Java reference outputs.

Regenerate fixture volumes with `make_phantom.py` from the repository root's
`pymacruise` directory. The ACE gray/white dimensions and analytic geometry
match `../ace/`; this directory separately commits all images needed by the
composed branch and the lesion/no-lesion segmentation variants.
