# LesionTOADS centroid-similarity formula reference

`java_reference.json` commits four deterministic two-channel/four-class
selected-method formula cases for the opening centroid-similarity section of
`LesionToadSegmentation.computeMemberships`. The cases exercise T1-specific
Sulcal-CSF weighting, the no-T1 fallback, lesion/black-hole/optimized-lesion
similarity exclusions, constructor-initialized black-hole similarity state,
and IEEE-754 results from zero standard deviations. Inputs and outputs are
stored as float32 bit patterns; NaNs use Java's canonical `Float.floatToIntBits`
encoding.

Regenerate the selected Java reference with
`scripts/generate_java_lesion_toads_similarity_reference.sh`. This validates
the coefficient setup only; no voxelwise membership or segmentation is
performed.
