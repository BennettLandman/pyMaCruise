# TOADS class-pair similarity prelude phantom

This small synthetic reference isolates the pairwise similarity-matrix
construction inside the selected VandyCRUISE
`MultipleToadSegmentation.computeExactMemberships` implementation. Four
classes (`OBJECT, LESION, LESION, OBJECT`) and two modalities exercise OO,
LL, and both mixed directions with deliberately asymmetric class and lesion
centroids and standard deviations. A second case sets two class deviations
to zero to retain Java IEEE NaN/infinity behavior.

`java_reference.json` was produced by
`scripts/generate_java_toads_pair_similarity_reference.sh` from the pinned
formula excerpt in `tests/java_reference/ToadsPairSimilarityReferenceHarness.java`.
The harness records float bits for the source coefficient, ordinary centroid
distance, and the mixed-pair lesion distance. For pairs where the source does
not assign `lesionDist`, the harness initializes it to NaN as a serialization
sentinel; this sentinel is not claimed to be a Java local-variable value.

This is a coefficient-prelude harness only. It does not validate the rest of
`computeExactMemberships`, membership updates, priors, or the TOADS pipeline.
