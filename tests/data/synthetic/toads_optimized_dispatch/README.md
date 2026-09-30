# Optimized TOADS dispatch fixture

This compact, deterministic fixture exercises one fixed-input membership pass
with the interleaved class sequence `OBJECT`, `OPTIMIZED_CSF`, `MASK`,
`OPTIMIZED_GM`, `OPTIMIZED_WM`. It has three intensity channels, a folded and
asymmetric ribbon mask, spatially varying intensity-field and PSF-bias maps,
and a true-RAS NIfTI affine with anisotropic millimeter spacing.

Regenerate the NIfTI inputs with:

```sh
python scripts/generate_toads_optimized_dispatch_fixture.py
```

The checked-in `java_reference.json` was produced by
`scripts/generate_java_toads_optimized_dispatch_reference.sh`. That script
verifies the selected VandyCRUISE sources by SHA-256 before compiling the
reference harness against the original `PowerTable.java`:

- `MultipleToadSegmentation.java`: `63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884`
- `PowerTable.java`: `5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678`

The reference contains exact float32 bit patterns for all output membership
voxels and the objective in four field/PSF combinations. The harness composes
only the explicit-input data terms and selected ordered membership-update
excerpt. It does not claim to implement full TOADS/MGDM processing or the
unhandled LESION, OUTLIER, OUTMASK, RAYLEIGH, auxiliary-map, optimized-update,
statistics-estimation, or preprocessing branches.
