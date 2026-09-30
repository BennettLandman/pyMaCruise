# LesionTOADS OUTLIER centroid phantom

`outlier_centroid_inputs.nii.gz` is a deterministic 5×4×3, 1-mm RAS phantom
with 12 float32 frames: hard labels, four class memberships, two image
channels, two correction fields, lesion/black-hole maps from the source
variance phantom, and an appended OUTLIER membership map. The map cycles
through 0, 0.25, 0.5, 0.75, and 1. The class-type vector selects OBJECT,
OUTLIER, OUTLIER, and OUTMASK; hard class 3 is intentionally absent.

Generate the image with `scripts/generate_lesion_outlier_centroids_fixture.py`
and its selected-method Java reference with
`scripts/generate_java_lesion_outlier_centroids_reference.sh`. The Java
excerpt and Python implementation compare raw float32 centroids for `none`,
`linked`, `damped`, `prior`, and fallback modes, with and without fields. This
validates OUTLIER membership suppression (`mem × (1 - outlierMembership)`) on
fixed hard labels only, not OUTLIER-map estimation or full LesionTOADS
iteration.
