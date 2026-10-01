# LesionTOADS initial loop-entry reference

This formula reference covers the `AlgorithmLesionToads.runAlgorithm()` state
transition immediately after the first standard-13 `computeMemberships()` call.
Its initial energy is pinned to the Java objective from
`lesion_standard13_volume/java_reference.json`; it models the subsequent
`setGMflag()`, `Eprev=energy`, `Eratio=1.0f`, iteration/stop initialization,
and the main-loop gate. No image payload is required for this scalar control
handoff; the upstream membership maps and objective are supplied by the
separate RAS synthetic-volume reference.

The selected source starts `propagateThinning()` when the gate is entered.
This seam stops immediately before that call. The existing general thinning
kernel has its own selected-source tests, but it has not been composed here
with the selected 13-class initialized state. Likewise, the subsequent
13-class membership call needs a separate `GMflag=true` source comparison;
the initial membership kernel intentionally models its earlier `false` state.
