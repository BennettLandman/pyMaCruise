# LesionTOADS segmentation configuration handoff

This fixture stores two-channel float32 intensity maxima and scales as exact
bit patterns. They are per-channel scalar statistics, not spatial images, so
JSON is the appropriate fixture format.

The selected wrapper branches to `setLesionWeight` or `reSetLesionWeight`,
then calls `setIntensityMax`, `setIntensityScale`, and `setGMflag` in that
order. The Java harness checks both lesion-weight branches, enabled GM flag,
array reference retention, and exact float bits. It does not validate how
preprocessing computes the intensity arrays.

`wrapper_configuration_source_excerpt.java.txt` preserves the selected
`AlgorithmLesionToads.runAlgorithm()` call sequence beginning with
`initCentroids`. The test verifies its pinned SHA-256 and checks that the
centroid initialization, lesion-weight branch, intensity references, and GM
flag calls appear in source order. The excerpt is independent reachability
evidence; the Java harness models the setters and float arrays.
