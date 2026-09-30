# Base-TOADS OUTMASK pass fixture

The committed synthetic inputs form an asymmetric, folded true-RAS 9×8×7
phantom with two channels and the ordered class sequence `MASK`, `OUTMASK`,
`OBJECT`. They include matching image/field/PSF-bias grids, active shapes,
all four closest-label slots, initial class memberships, and an initial
outlier-membership map. Every class exercises DIRECT, NEIGHBOR, and PRIOR
relations; the pass also uses nonzero smoothing and a nonzero atlas prior.

Regenerate the NIfTI files with:

```sh
python scripts/generate_toads_outmask_fixture.py
```

Regenerate the pinned Java reference with:

```sh
scripts/generate_java_toads_outmask_reference.sh
```

The harness pins `MultipleToadSegmentation.java` SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`. It
matches class data energies, OUTMASK auxiliary energy, updated class and
outlier memberships, and objective bits in field/PSF on/off combinations.
Its unbraced Java LESION smoothing conditional reproduces the selected
OUTMASK auxiliary-smoothing skip; the braced prior term is still applied to
the auxiliary energy. This fixture validates one fixed-input pass, not full
TOADS/MGDM or other lesion/outlier families.
