# LesionTOADS post-thinning registration fixture

This fixture exercises the selected VandyCRUISE `linear_fast_marching` source
order through `computeMemberships(GMflag=true)` → energy-ratio calculation →
`registerShapes()`. It stops at registration. It does not model or assert the
subsequent convergence decision or any later loop branch.

`registration_inputs_ras.nii.gz` is a synthetic 17 × 17 × 17 RAS-grid volume
with 26 float32 frames: the 13 post-thinning class memberships followed by
their 13 class shape-prior fields. The membership frames are from the folded
synthetic P39ak fixture and the fixed registration priors are supplied as
fixture inputs. No actual imaging data is used. Registration runs at 1 mm
isotropic image/template resolutions from the fixed six-parameter rigid
transform recorded in `fixture.json`, with `subsample=2` and
`init_align_iter=10`.

The selected membership update has a NaN objective/energy ratio caused by
source-faithful zero-tensor divisions, while its class membership arrays are
finite. The VandyCRUISE wrapper calls `registerShapes()` before its later
convergence check, so the NaN ratio does not suppress this registration call.
This fixture makes no claim about that later stop decision.

Regenerate the NIfTI fixture and Java reference with:

```sh
sh scripts/generate_lesion_post_thinning_registration_reference.sh
```

`java_reference.json` records the selected Java registration optimizer trace,
final transform and energy, rotation, transformed support bounds, source-file
hashes, and harness hash. The harness uses the checked-in Java numerics,
transform, interpolation, and rotation source snapshots under
`tests/java_reference/upstream/`.
