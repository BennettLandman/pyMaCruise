# Default TOADS membership-to-registration handoff phantom

This handoff uses the folded, asymmetric 9×8×7 two-channel true-RAS subject
inputs in `../toads_default_atlas_dispatch/` and a dedicated 11×10×9
`atlas_shapes.nii.gz` template prior. The larger template grid is intentional:
the source keeps subject membership dimensions separate from atlas dimensions.
The fixture is synthetic; it contains no subject imaging data. The membership pass uses the selected default atlas
order (MASK, 2 optimized CSF, 5 optimized GM, 3 optimized WM), then its newly
updated maps are passed to the rigid atlas-registration adapter when
`align_shapes=True`.

The Java call-sequence excerpt records the selected wrapper order:
membership and energy ratio, optional shape registration, convergence and
iteration-limit checks, then growth and the stop/continue state updates.
`AlgorithmMultipleToads` defaults `alignShapes` to true and initializes this
registration mode as rigid SINGLE-transform registration. `registerShapes`
resets lambda, turn count, and positive-update count; repeats the gradient
descent while `n < itMax`, `diff > minEdiff`, and `lambda > minLambda`; then
refreshes rotation and transformed shape bounds. The Python adapter composes
the first two numerical stages; existing P18cc/P18ck selected-Java checks
cover the registration optimizer separately. This is a fixed-input handoff,
not a full-wrapper Java execution or segmentation-equivalence claim.

Generate the dedicated shape map with:

```sh
python scripts/generate_toads_default_registration_handoff_fixture.py
```

Generate the wrapper-order reference with:

```sh
scripts/generate_java_toads_default_registration_handoff_reference.sh
```
