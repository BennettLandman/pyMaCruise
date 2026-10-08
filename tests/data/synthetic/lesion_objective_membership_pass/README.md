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
records bitwise results for the supported membership maps and shared
denominator. The test checks those maps, the input-data digest, and that all
six boundary faces are unchanged. The harness energy/objective arrays are not
source-equivalence evidence: the harness resets per-voxel energy, while the
selected Java method uses method-local accumulators across its full x/y/z scan.
The source-faithful Python option and zero-shape regression now preserve that
behavior, but a direct Java composite reference for its objective remains open.

All spatial and data-energy arrays are fixed inputs; the harness compares the
supported interior membership, smoothing, prior, and normalization maps.
Shape-prior values, similarities, relations, smoothing coefficients, and
class order are fixed in the test/harness. All classes stay active in this
fixture. A separate source-faithful regression covers the source's stale
`energy[k]` value when shape suppresses a class. Unsupported types, input
energy construction, shape/topology/distance generation, and full iteration
or workflow behavior are outside this fixture.

Regenerate the RAS NIfTI and selected-Java reference with:

```sh
.venv/bin/python tests/data/synthetic/lesion_objective_membership_pass/generate_fixture.py
sh scripts/generate_java_lesion_objective_membership_pass_reference.sh
```
