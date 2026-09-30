# TOADS dura-removal output reference

This asymmetric synthetic reference exercises the selected
`MultipleToadSegmentation.generateDuraRemovalOutputs()` output arrays, the
four scalar-image uncrop/import branch in `AlgorithmMultipleToads`, and the
selected `MedicAlgorithmToads08` wrapper's unconditional multiplication by 255. The
selected-method Java excerpt harness was compiled and run to create
`java_reference.json`; it does not instantiate the surrounding MIPAV/JIST
pipeline or serialize NIfTI images.

## Source provenance

Selected checkout: `/Users/landmaba/dev/mipav/toads-cruise-scm-2026-09-26/VandyCRUISE/`

- `src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`, selected
  method `generateDuraRemovalOutputs`, lines 3245–3299; full-file SHA-256:
  `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
- `src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java`, selected
  `Cruise_input` / `DuraRemoval_input` branch, lines 884–924; full-file
  SHA-256:
  `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`.
- `src/edu/jhmi/rad/medic/methods/MultipleToadPreprocess.java`, selected
  `uncropAndBuffer(float[][][], float)` and `uncropAndBuffer(byte[][][], byte)`
  overloads, lines 478–517; full-file
  SHA-256:
  `31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7`.
- `toads-cruise/src/edu/jhu/ece/iacl/plugins/classification/MedicAlgorithmToads08.java`,
  selected output branch lines 365–385; full-file
  SHA-256:
  `1153b5aa063a3d3f97904494c5c5097c8546c657810d1f3e8044562d24518261`.
- Harness SHA-256:
  `d2af265c37a632fafa99c9061851c2eb4ffa9f930054a6b5d43624d434039b7d`.
- Generated reference SHA-256:
  `85dc342d530d82c07edbc2bfabf6c5243ceedd7472717a99d2dca62bb4a73437`.

Regenerate the reference with:

```sh
javac -d /tmp/toads-dura-output-reference tests/java_reference/ToadsDuraRemovalOutputsReferenceHarness.java
java -cp /tmp/toads-dura-output-reference ToadsDuraRemovalOutputsReferenceHarness
```

## Fixture and preserved behavior

The crop grid is `[5,4,4]`, restored into an original `[4,3,3]` extent using
bounds `[-1,1,0,1,-1,0]` and one-voxel padding. The eight classes include both
CSF aliases, cortical GM, both recognized cerebellar spellings, the literal
`CerebellarGM)` typo, a white-matter class, and background. Memberships and
labels vary by voxel so class order, sums, label inclusion, clipping, and
x-fastest output buffer indexing are observable.

The selected source quirks are deliberate in both harness and Python:

- CSF aliases are assigned in class order, so the later matching alias wins.
- GM soft membership includes correctly spelled `CerebellarGM`.
- WM fill excludes `CerebellarGM`, but not the literal `CerebellarGM)` typo.
- The hard mask sets Cerebellum-GM and `CerebellarGM)` labels in its GM loop;
  the correctly spelled `CerebellarGM` label does not pass that comparison.
- The algorithm adapter creates four scalar FLOAT images and restores each
  product with zero background. This is the raw buffer stage.
- The hard classification is restored separately into a byte buffer using
  background value 1.
- The selected lower-case-option `MedicAlgorithmToads08` source comments out
  both the `scaleData` parameter and its guard, then unconditionally multiplies
  CSF, GM, WM fill and WM mask by 255 after image import. Python returns the
  raw buffers and the always-scaled wrapper buffers as separate stages.

The harness calls a small float-array multiplication helper to model the
`ImageDataMath.scaleFloatValue(...,255.0f)` operation named by the wrapper. The
ImageDataMath dependency implementation is not present in the available
source checkout, so this comparison validates the selected wrapper's call
sites/default and IEEE float32 multiply behavior, not the external helper's
implementation or MIPAV's image object lifecycle.

`tests/test_toads_dura_outputs.py` pins the harness and generated-reference
hashes and compares raw buffers, byte classification, and scaled buffers
exactly. The fixture validates these selected array computations and packing
only, not MIPAV `ModelImage` import, NIfTI headers/serialization, or a complete
TOADS execution.
