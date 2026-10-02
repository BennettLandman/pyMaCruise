# Carried-degree exact TOADS iteration (P18cyba)

## Ported behavior

The selected next-loop continuation now carries the polynomial degree from
the prior iteration through field fitting, centroid and variance refresh,
image-variance and centroid-similarity recomputation, and the next exact
membership pass. A qualifying energy ratio advances the degree by one, up to
degree three. Degree zero preserves the attached field and proceeds with the
existing field. A positive degree refits every iteration, including degree
three when the degree no longer increments; this follows the wrapper's
unconditional `computeCorrectionField(deg)` call after its degree gate.

The single-channel 11-class profile, fixed registration/shape/relation state,
and membership parameters remain explicit inputs. This does not assemble the
outer TOADS workflow or estimate those upstream inputs.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_post_thinning_exact_field_ramp.py tests/test_toads_post_thinning_exact_field.py tests/test_toads_post_thinning_exact_iteration.py tests/test_toads_p18cza_carried_degree_iteration.py
```

Result: **13 passed, 0 failed, 0 skipped in 178.40 seconds**. The tested tree
was based on `fc30d686a3aaae92036071ce71ff801b27eb5989`; its implementation
and test changes are included in the commit containing this report. The only
source edit after the run corrected a docstring to match the tested behavior.

The four carried-degree cases cover 0→1, 1→2, 2→3, and a degree-3 plateau.
The degree-1/2/3 fields match the selected Java formula harness within
`2e-7`; the Java formula iteration harness then checks centroids, deviations,
similarity, image variance, membership maps, denominator, and objective using
the refitted field. The plateau case confirms degree 3 is still refitted.
Existing ramp tests also check the degree-zero no-op and its high-ratio
control. The selected 21³ folded-cortex input is RAS and committed as NIfTI.

## Scope limits

The composition supports one channel in this adapter. Field and statistics
kernels have separate multichannel implementations/tests. Registration,
transformed-shape refresh, topology relation generation, growing, and outer
loop termination remain separate stages. This is selected synthetic
method/composition evidence, not whole-TOADS equivalence.
