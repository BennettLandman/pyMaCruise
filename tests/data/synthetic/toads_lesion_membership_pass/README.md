# Base TOADS OBJECT+LESION synthetic membership fixture

`inputs_ras.nii.gz` is a deterministic 5×4×5×16 float32 RAS+ NIfTI fixture. It
contains two image channels; fixed multiplicative fields and PSF biases; initial
OBJECT/LESION class memberships; a separate lesion-membership map; positive
shape values; a categorical segmentation map; and the four selected closest
class labels. It is asymmetric and folded by modular spatial patterns, with
dedicated zero-residual voxels for reciprocal-threshold cases. It contains no
patient data.

The Java reference harness pins
`VandyCRUISE MultipleToadSegmentation.computeExactMemberships()` at SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`. It covers
OBJECT+LESION class orders, segmentation-label matches and mismatches, field
and PSF-bias modes, independent smoothing/prior switches, direct/neighbor/prior
relations, source direct-relation precedence, raw terms, updated maps, the
normalization denominator, and objective accumulation.

Regenerate the fixture and Java excerpt output from the repository root with:

```sh
sh scripts/generate_java_toads_lesion_membership_reference.sh
```

This validates one fixed-input exact-membership branch. It does not validate
statistics estimation, PSF-bias estimation, topology checks, inactive-shape
stale-energy behavior, optimized lesion types, distance relations, or full
TOADS orchestration.
