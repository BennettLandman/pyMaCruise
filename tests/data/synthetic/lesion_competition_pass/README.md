# LesionTOADS competition-pass phantom

`competition_inputs.nii.gz` is a 17×17×17, 4-D float32 synthetic RAS image
with frames for the categorical source labels, four class memberships,
lesion membership, black-hole membership, and prior ordering. It contains a
folded asymmetric ribbon with a narrow branch and internal lesion/black-hole
labels. Membership values are deliberately sparse in the background so the
selected Java queue reaches an empty state without scheduling self-transitions
at the untested volume border.

`java_reference.json` was emitted by the selected Java excerpts for
`initCompetitionLabels()` and `propagateCompetition()` over the serialized
NIfTI input. The harness uses the selected Java `BinaryHeap4D` and direct
`LesionToadSegmentation.isHomeomorphicSegmentation()` logic with the bundled
18/6 critical-point LUT. It records initial queue entries, every pop and
requeue, stale and topology outcomes, and hashes of the complete final state.
It runs both lesion object-type variants (`LESION` and `BLACKHOLE`).

`java_topology_reference.json` independently pins 140 generated local
neighborhood decisions from the selected LesionTOADS topology method and LUT;
the Python topology helper must match those decisions before it is reused in
the pass.

This is one source-scoped competition queue pass from explicit fixed state.
It does not validate initialization/preprocessing, membership/statistics
estimation, atlas alignment, the wrapper loop, growing, or complete
LesionTOADS behavior. Regenerate the NIfTI with
`scripts/generate_lesion_competition_fixture.py`; regenerate references by
compiling the harnesses with the pinned Java heap/LUT sources and running the
competition harness on that exact `.nii.gz` file.
