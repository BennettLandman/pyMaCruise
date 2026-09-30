# LesionTOADS objective-reporting membership pass fixture

`inputs_ras.nii.gz` is a deterministic 5×4×5 float32 RAS NIfTI with an
asymmetric folded value pattern. Its 18 channels are grouped as follows:

1. Initial class memberships: OBJECT, LESION, BLACKHOLE, OUTMASK.
2. Initial lesion, black-hole, and outlier memberships.
3. Raw class data terms: OBJECT, LESION, BLACKHOLE, OUTMASK.
4. Raw lesion, black-hole, and OUTMASK data terms.
5. Positive transformed-shape values for all four classes.

The Java excerpt harness pins selected
`LesionToadSegmentation.computeMemberships()` source SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99` and
records bitwise results for every returned map, the shared denominator, class
and auxiliary objective-energy maps, the running objective, and final total.
The test also checks the input-data digest and that all six boundary faces are
unchanged.

All spatial and data-energy arrays are fixed inputs; the pass composes only
the supported interior membership, smoothing, prior, normalization, and
objective behavior. Shape-prior values, similarities, relations, smoothing
coefficients, and class order are fixed in the test/harness. All classes stay
active in the fixture so it does not imply support for the source's stale
`energy[k]` value when shape suppresses a class. Unsupported types, input
energy construction, shape/topology/distance generation, and full iteration
or workflow behavior are outside this fixture.

Regenerate the RAS NIfTI and selected-Java reference with:

```sh
.venv/bin/python tests/data/synthetic/lesion_objective_membership_pass/generate_fixture.py
sh scripts/generate_java_lesion_objective_membership_pass_reference.sh
```
