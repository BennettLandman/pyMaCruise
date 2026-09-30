# LesionTOADS topology-changing thinning phantom

`topology_thinning_inputs.nii.gz` is a deterministic 11×11×11, 1-mm RAS
synthetic volume packed with 14 frames: starting segmentation, skeleton, three
class memberships, lesion and black-hole maps, transformed shape maps, three
availability maps, and an ordering field. It uses a folded, mostly coherent
three-class pattern with 5% deterministic local label perturbations. Per-voxel
memberships favor the cyclic next class (0→1, 1→2, 2→0), so the thinning pass
must make source-gated class changes rather than only remove boundary points.

The Java and Python generators use the same documented 31-bit linear
congruential sequence (`state = (1103515245*state + 12345) & 0x7fffffff`),
seed 1, x/y/z traversal, and 3-voxel coarse folds. Lesion/black-hole maps and
transformed shape priors are zero, allowing the direct six-face target path to
be observed without conflating shape fallback. The minimum-volume ratio is
0.4 with the selected source floor of 100; the queue limit is 10.

The selected-method Java harness calls the bundled source `CriticalPointLUT`
through the already source-compared 18/6 topology predicate and uses the copied
`BinaryHeap4D` and `Numerics.max`. It produces a 706-event trace: 9 accepted and
69 rejected topology decisions, 9 final changed labels, and 5 newly inserted
neighbors. Python matches every popped event, candidate and topology decision,
neighbor insertion, final array, volume float32 bit, and remaining queue.

Regenerate with `scripts/generate_lesion_toads_topology_thinning_fixture.py`
and `scripts/generate_java_lesion_toads_topology_thinning_reference.sh`.
This validates one fixed-input selected thinning pass, not LesionTOADS atlas
sampling, membership estimation, registration, or complete segmentation.
