# LesionTOADS variance phantom

`variance_inputs.nii.gz` is a deterministic 5×4×3, 1-mm RAS synthetic volume
with 11 float32 frames: hard labels, four class memberships, two image
channels, two correction fields, and lesion/black-hole membership maps.
The hard labels contain ordinary, LESION, and BLACKHOLE classes plus one
absent class. The maps exercise lesion and black-hole winners, exact ties, and
locations where both auxiliary memberships are zero.

The selected-method Java excerpt and Python implementation match raw float32
bits for class, lesion, and black-hole standard deviations with and without
field correction. The Java reference deliberately retains the source's
BLACKHOLE winner behavior: it weights the black-hole residual by the lesion
membership value. This source behavior is documented and reproduced, not
corrected silently. Class standard deviation for the absent class and
auxiliary standard deviations with no selected support are zero.

Regenerate the NIfTI with `scripts/generate_lesion_variance_fixture.py` and
the reference with `scripts/generate_java_lesion_variances_reference.sh`.
This validates the variance updater for supplied hard labels and membership
maps; it does not validate centroid-category assignment or full membership
iteration.
