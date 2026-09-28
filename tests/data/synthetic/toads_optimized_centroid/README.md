# Optimized TOADS centroid phantom

This source comparison reuses the committed 7×6×5×3 RAS `.nii.gz` images and
intensity field in `toads_optimized_likelihood/`, and the first three class
membership channels in `toads_direct_membership/`. No duplicate image files
are stored. Initial CSF/GM/WM centroids and exponents are recorded in
`java_reference.json`; one class/modality centroid exactly equals a synthetic
image value to exercise the source zero-distance factor path.

The Java harness compiles against the preserved actual `PowerTable.java`
class. It implements the selected optimized-CSF/GM/WM weighted mean and the
MaCRUISE wrapper's damped centroid mode (`centroidSmoothness=0.1f`) with the
intensity field disabled and enabled. Source hashes are stored in the JSON.
It does not cover object, lesion, outlier, Rayleigh, linked/prior centroid
modes, or the full segmentation loop.

Regenerate the Java reference with
`scripts/generate_java_toads_optimized_centroid_reference.sh`.
