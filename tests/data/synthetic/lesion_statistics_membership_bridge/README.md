# P39zzzzs LesionTOADS statistics-to-membership bridge

`bridge_inputs_ras.nii.gz` is a deterministic 6×5×4, two-channel float32
RAS image containing hard labels, current class memberships, lesion,
black-hole and outlier maps, correction fields, and positive class shapes.
Its formulas make folded, asymmetric boundaries and distinct data values on
all axes. The same input is evaluated both with correction fields disabled
and enabled.

`java_reference.json` records raw float32 bits from
`LesionStatisticsMembershipBridgeReferenceHarness.java`. The harness is a
selected-source formula composition of wrapper order:

1. `computeCentroids()`, `computeVariances()`, `computeChannelWeights()`
2. Recomputed centroid similarities and raw OBJECT/LESION/BLACKHOLE energy
   maps, plus the selected OUTMASK `maskDistance` and outlier-distance terms
3. The `computeMemberships()` map/objective pass
4. Wrapper `Eratio = abs(2*(Eprev-E)/(Eprev+E))`

The selected `LesionToadSegmentation.java` SHA-256 is
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`; the
`AlgorithmLesionToads.java` wrapper SHA-256 is
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`.
The reference pins the harness hash as well. Rebuild the NIfTI and Java JSON
with `scripts/generate_java_lesion_statistics_membership_bridge_reference.sh`.

The bridge accepts fixed post-thinning labels, memberships, fields, shapes,
priors, relations, and previous energy. It does not estimate correction
fields, thin or change topology, register atlas shapes, decide convergence,
grow labels, or execute a complete LesionTOADS iteration. The fixture covers
OBJECT, LESION, BLACKHOLE, and OUTMASK in this order and direct/prior
relations; other object-type families and source distance-modulation flags
are excluded.
