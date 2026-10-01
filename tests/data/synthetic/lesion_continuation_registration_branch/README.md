# P39zzzzw registration/branch composition fixture

This bounded composition reuses the synthetic asymmetric true-RAS 15×15×15
multi-channel input at `../lesion_continuing_iteration/iteration_inputs_ras.nii.gz`.
That fixture drives the P39v continuing iteration and has its own generator,
provenance, and method-level Java comparisons. No clinical images are used.

The Java wrapper excerpt fixes the P39v profile's resulting previous and
membership objectives at float32 values 49.18265151977539 and
26.26949119567871. The selected wrapper calculates its float32 energy ratio,
registers shapes, checks convergence before the absolute signed iteration
ceiling, grows, then refreshes relations when continuing.

The optimized registration kernels and post-registration growth/relation
kernels keep their own selected-Java method harness evidence in
`test_lesion_registration_handoff.py` and
`test_lesion_post_registration_branch.py`. This fixture pins their new
composition boundary and wrapper controls; it does not present a new Java
optimizer run on the P39v-produced membership maps as having occurred.
