# Optimized-tissue plus OUTMASK fixture

This committed synthetic input set is a 9×8×7 asymmetric, folded, two-channel
phantom in true RAS orientation. Class order is `OPTIMIZED_GM`, `OUTMASK`,
`MASK`, `OBJECT`, `OPTIMIZED_CSF`, `OPTIMIZED_WM`. The six classes use active
folded shape maps and labels exercising direct, neighbor, and prior routes.
The image, field, and PSF bias share the same anisotropic affine; initial
class and outlier memberships are separate inputs.

Regenerate the NIfTI fixture with:

```sh
python scripts/generate_toads_optimized_outmask_fixture.py
```

Regenerate the selected-Java reference with:

```sh
scripts/generate_java_toads_optimized_outmask_reference.sh
```

The harness pins `MultipleToadSegmentation.java` SHA-256
`63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884` and
`PowerTable.java` SHA-256
`5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678`. Its
fixed-input excerpt compares all class and auxiliary energy, memberships,
outlier map, denominator, and objective float32 bits in field/PSF and
smoothing/prior on/off combinations. It includes the selected OUTMASK
dangling-else smoothing behavior and category-specific optimized exponents;
it does not claim full TOADS/MGDM workflow equivalence.
