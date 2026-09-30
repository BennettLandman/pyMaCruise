# LesionTOADS thinning-frontier phantom

The six committed `.nii.gz` inputs form an asymmetric 4×4×4, 1-mm RAS
phantom with three classes (`OBJECT`, `LESION`, `BLACKHOLE`). Array values are
deterministic integer fractions of 1/32 or 1/64; class labels are
`(x+y+z) mod 3`, with the six face neighbors of `(1,1,1)` reassigned to class
0 so that this interior voxel is not a boundary candidate. The initial
skeleton contains sentinel 7 on the outer planes and -1 in the interior.
Initial availability alternates by coordinate and class, making it possible
to verify that only the interior state is replaced.

`java_reference.json` records both `useLesionMax` modes from the selected
`LesionToadSegmentation.initThinningLabels` source-method excerpt and queue
removal order from the copied selected-source `BinaryHeap4D`. The harness
asserts the segmentation source SHA-256; the copied heap SHA-256 is
`9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a`.
The Java heap is used because equal priorities follow its insertion/removal
path rather than a stable sorting rule.

This fixture checks reset class volumes, six-face boundary discovery,
interior skeleton/availability updates, source-order queue insertions, and
the complete max-heap pop sequence. It does not exercise topology decisions
or removal events in `propagateThinning`.
