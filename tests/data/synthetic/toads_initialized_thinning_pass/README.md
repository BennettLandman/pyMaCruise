# Initialized TOADS thinning pass

This 12×12×12 four-class analytic phantom runs selected VandyCRUISE
initThinningLabels and then propagateThinning to queue exhaustion. Labels
are (x + 2y + 3z) % 4; memberships use
(((5x + 7y + 11z + 3k) % 29) + 1) / 30 with float32 arithmetic. It contains
no subject images.

The four initialized class volumes are each above 100. With
volumeRatio=0.8, class 2 uses the 100-voxel floor and the other classes use
volumeRatio * volume[k]; all minima are recorded as float32 bit patterns.
The pass uses limit=3, an all-false relation matrix, and zero shape priors so
this fixture focuses on initialization and source-ordered direct-candidate
thinning. The Java reference and Python output match the initialization
volumes and boundary counts, all 1,000 popped heap payloads/priorities, final
segmentation/skeleton/availability, post-pass class volumes, and empty
residual queue exactly.

The selected source's boundary counter includes all non-face neighbors, so
its endpoint branch cannot occur on a 3×3×3 neighborhood; the port preserves
that documented behavior. Regenerate with
scripts/generate_java_toads_initialized_thinning_pass_reference.sh. This
validates one initialized synthetic thinning pass, not atlas registration,
growing, Dura, or complete TOADS segmentation.
