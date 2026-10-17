# AtlasEM direct Demon pyramid optimizer comparison — 2026-10-17

## Finding

The initial AtlasEM pyramid registration now matches the selected Java
`DemonToadDeformableAtlas.registerShapesPyramid()` optimizer on a controlled
synthetic input. Both Java and Python execute the same pyramid levels and
produce bit-identical values for all six final float32 rigid-transform
parameters at one- and two-iteration caps.

This comparison exercises the selected Java class directly. It does not
substitute the related TOADS optimizer or infer whole-workflow equivalence.

## Synthetic case and reference

- Generated anisotropic three-dimensional image, membership, and class-prior
  arrays on a `16 × 16 × 16` grid.
- Three pyramid levels are requested; the selected Java method optimizes at
  scales 4 and 2, with no additional scale-1 pass.
- The harness calls the selected class's `initShapeRegistration(..., levels)`
  and `registerShapesPyramid()` methods. Python uses the AtlasEM pyramid
  adapter with the same arrays and iteration cap.
- The comparison checks the exact float32 bit patterns of all six rigid
  transform parameters. The test exercises caps of one and two updates.

| Iteration cap | Transform comparison | Result |
|---:|---|---|
| 1 | Six Java/Python float32 parameter bit patterns | Bit-identical |
| 2 | Six Java/Python float32 parameter bit patterns | Bit-identical |

This is narrow synthetic optimizer parity. It does not validate descriptor
loading, the unavailable selected `ibsr8obj` payload, the subsequent
post-recrop stage in context, AtlasEM convergence on human data, or the
AtlasEM-to-LesionTOADS PLUS route.

## Validation

Focused direct comparison: **4 passed in 4.10 s**.

The full suite ran on code/test commit
`f179290210fb2f0c9cdb193833631dd473081799`: **1,650 passed, 4 skipped,
0 failed, 21 expected orientation warnings in 719.57 s**. Peak sampled RSS
was **14.822 GB** under the 16 GB watchdog; the watchdog did not terminate the
run. See the [full-suite report](full-pytest-2026-10-17-atlasem-direct-demon-pyramid.md)
for the exact command.
