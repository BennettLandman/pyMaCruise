# TOADS closest-label propagation fixture

This asymmetric 9×8×7 synthetic six-class segmentation exercises initial
face-neighbor ordering and synchronous diffusion of the selected four-slot
`closestLabel` map. Java and Python references are generated for limits 0, 1,
and 3. Every internal class/slot is compared, including `EMPTY=-1` and the
unchanged volume border. There are no medical images in this fixture.

Regenerate with `scripts/generate_java_toads_closest_labels_reference.sh`.
The harness is a selected-method excerpt of
`MultipleToadSegmentation.computeClosestLabels`; it is not full TOADS
segmentation equivalence. Source file SHA-256:
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
