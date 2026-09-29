# TOADS local thinning transition fixture

This deterministic fixture validates one interior queue item from the selected
VandyCRUISE `MultipleToadSegmentation.propagateThinning` method. The inputs are
small label and skeleton neighborhoods, four-class memberships/relation flags,
transformed shape priors, and source class-volume thresholds. There are no
image payloads in this formula-level fixture.

Nine cases exercise unchanged-best class, exact volume threshold, a stale
skeleton item, no connected replacement, direct topology acceptance and
rejection, relation-connected replacement acceptance and rejection, and a
shape-prior tie resolved in source face-neighbor order. Discrete labels,
skeleton, candidate, topology outcome, availability flags, and float32 volume
updates are recorded from a Java formula harness that invokes the selected
Java 18/6 topology-reference method and committed lookup table.

Regenerate the expected arrays with
`scripts/generate_java_toads_thinning_transition_reference.sh`. The harness
uses the upstream `CriticalPointLUT` implementation and the bundled
`critical186LUT.raw.gz`; neither is a runtime JVM dependency. This is a local
transition only. It does not validate boundary-queue priority/tie behavior,
initialization over a volume, repeated thinning, subsequent centroid or
membership updates, or full TOADS segmentation.
