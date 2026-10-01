# P18cye synthetic reference fixture

This fixture isolates `computeApproxMemberships()` through the alternate
`AlgorithmMultipleToads` constructor path that accepts `approx=4`, at one fixed
TOADS state. The selected `MedicAlgorithmToads08` wrapper calls the shorter
constructor, which explicitly sets `useApprox=false`; this fixture does not
claim that the wrapper runs this approximate path. It is synthetic only.
Images, fields, sparse
candidate memberships, and best-class slots start from the earlier synthetic
P18cyc/P18cyd state, reduced to one channel; folded transformed shapes and relation membership
are added deterministically. All NIfTI inputs have a true-RAS affine with
0.85, 1.15, and 1.55 mm voxel spacing.

The Java reference harness pins the method, alternate constructor, default
descriptor, `Numerics`, and `PowerTable` source hashes:

- `MultipleToadSegmentation.java`: `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`
- `AlgorithmMultipleToads.java`: `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`
- `MedicAlgorithmToads08.java`: `d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115`
- selected atlas descriptor: `d43efa5a002b0c1a838449ba7ca1d3bf35fab5e3b575474a71aadfe475dcc257`
- `Numerics.java`: `13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb`
- `PowerTable.java`: `5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678`

Regenerate the inputs and source-formula reference with:

```sh
.venv/bin/python tests/data/synthetic/toads_default_approximate_pass/generate_fixture.py
sh tests/java_reference/generate_toads_default_approximate_membership_pass_reference.sh
```

The harness compares image-variance and centroid-similarity prelude outputs,
four updated memberships and class IDs, per-voxel denominators, and the total
objective as exact float32 bits. It runs both fielded and raw intensity paths.
The fixture contains exact candidate-score ties, zero and near-zero costs,
voxels where every class is prior-ineligible, and a first-swept interior voxel
whose six face neighbors can be perturbed independently. The untargeted image
faces also verify that the method leaves border membership state unchanged.

This validates only the fixed-input default 11-class, approx=4 alternate
constructor method path. It supplies centroids, class deviations/weights, transformed atlas
shapes, and relation state. It does not perform atlas initialization,
inhomogeneity fitting, topology/thinning, registration, or outer-loop control.
