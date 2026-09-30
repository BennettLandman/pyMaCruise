# TOADS Rayleigh distance fixture

This fixture stores 17 `(val, sig)` float32 pairs in `rayleigh_inputs.nii.gz`
with a true-RAS affine and anisotropic millimeter spacing. The pairs cover
positive and negative zero, negative and ordinary positive values, zero and
negative sigma, subnormals, large finite values, NaN, and positive/negative
infinity. NIfTI dimension four contains the two channels.

Regenerate the inputs with:

```sh
python scripts/generate_toads_rayleigh_fixture.py
```

Regenerate `java_reference.json` with:

```sh
scripts/generate_java_toads_rayleigh_distance_reference.sh
```

The harness reproduces the selected `rayleighDistance(float,float)` method
from `MultipleToadSegmentation.java`, SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
Reference outputs use Java `Float.floatToIntBits`, including canonical NaN
bits, and preserve the source's `val <= 0` return of exactly `1.0f`.
