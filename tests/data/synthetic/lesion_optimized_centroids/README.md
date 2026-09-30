# LesionTOADS optimized tissue centroid phantom

`optimized_centroid_inputs.nii.gz` is a deterministic 5×4×3, 1-mm RAS
phantom with 12 float32 frames: labels, five class memberships, two images,
two fields, and the two unused auxiliary maps retained from the base
variance phantom. Its classes are OBJECT, optimized CSF, optimized GM,
optimized WM, and an absent OBJECT class. The channels are T1_MPRAGE and
FLAIR, so optimized CSF is PowerTable-weighted only in T1_MPRAGE as in the
selected source branch.

The fixture includes an exact zero residual (with a unit field) and a small
uncorrected residual whose square is below the selected `ZERO=1e-10` cutoff.
Other optimized residuals exercise cached PowerTable weights and the direct
out-of-range behavior. Generate the NIfTI with
`scripts/generate_lesion_optimized_centroids_fixture.py` and the selected
Java reference with
`scripts/generate_java_lesion_optimized_centroids_reference.sh`. The harness
compiles against the preserved upstream `PowerTable.java` and compares raw
float32 centroids for all supported update modes, with and without fields.
This validates optimized CSF/GM/WM centroid weights on fixed maps only;
optimized LESION/BLACKHOLE and the complete membership iteration remain
separate.
