# LesionTOADS optimized auxiliary-centroid phantom

`optimized_auxiliary_centroid_inputs.nii.gz` is a derivative of the committed
5×4×3, 1-mm RAS LESION/BLACKHOLE phantom. It retains the five class-membership
frames, two images, two fields, lesion map, and black-hole map, with classes
1 and 2 typed as OPTIMIZED_LESION and OPTIMIZED_BLACKHOLE. Added cases include
a tiny winning class membership whose squared-membership intensity is below
`ZERO`, and an optimized lesion winner whose uncorrected lookup argument is
above the PowerTable cache maximum.

Generate the image with
`scripts/generate_lesion_optimized_auxiliary_centroids_fixture.py` and the
Java reference with
`scripts/generate_java_lesion_optimized_auxiliary_centroids_reference.sh`.
The reference compiles against the preserved upstream `PowerTable.java` and
compares class, lesion, and black-hole centroid float32 bits across all five
centroid modes, with and without field correction. The unusual source factor
is computed from the winning membership squared times corrected intensity,
then shared by the winning branch's class/auxiliary accumulator. This
validates fixed-map centroid routing only.
