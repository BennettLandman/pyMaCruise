# Initial LesionTOADS membership to first thinning pass

This folded, asymmetric 17×17×17, two-channel synthetic fixture is stored as
`initial_thinning_inputs_ras.nii.gz` on a 1 mm RAS+ affine. Its 64 frames are
named in `fixture.json`: T1 and FLAIR inputs, 13 transformed-shape maps, 13
direct-relation flags, the supplied class initialization, initial skeleton,
availability, ordering, three fixed distance maps, and zero-seeded membership
and auxiliary maps. The class initialization follows a folded interface
`z = 8 + ((3*x + 5*y) mod 5)` between classes 1 and 12.

The test runs initial standard-13 `computeMemberships()` with `GMflag=false`,
passes that result through the selected loop-entry transition (`setGMflag`,
initial energy and ratio, iteration-one gate), and immediately feeds its maps
to `initThinningLabels()` and `propagateThinning()`. The remaining statistics,
shape maps, class initialization, and atlas state are fixed synthetic inputs.
The profile marks class 10 as LESION, uses `volumeRatio=0.4`, a 1000 priority
limit, the selected 18/6 topology LUT, and the Java 4D max heap.

The pinned Java method harness compares every initial membership, energy,
relation, denominator, objective and running-objective map; initial and minimum
class volumes; all inserted queue entries and pop/update events; topology
decisions; remaining queue; and final segmentation, skeleton, availability,
ordering, and class volumes. Its FLAIR lesion term calls the supplied
`PowerTable` implementation using the selected cached float exponent `1.2f`
(bits `0x3f99999a`); it does not use precomputed fixture-specific power values.
The fixture exceeds the 100-voxel minimum for
classes 1 and 12 and exercises both accepted and rejected topology updates.
This ends after the first thinning pass. It does not run the later
`GMflag=true` membership update, centroid/variance refresh, registration,
growth, or a complete LesionTOADS workflow.

Regenerate the input with
`scripts/generate_lesion_initial_thinning_composition_fixture.py` and the Java
reference with
`scripts/generate_java_lesion_initial_thinning_composition_reference.sh`.
