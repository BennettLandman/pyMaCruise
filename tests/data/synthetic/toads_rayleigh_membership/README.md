# TOADS homogeneous RAYLEIGH membership-pass fixture

This asymmetric, folded 9×8×7 phantom contains three-channel image, field, and
PSF-bias volumes; a folded ribbon; explicit initial memberships and active
shape values; and four-slot closest-label relations. Every volume has the
same true-RAS affine and anisotropic millimeter spacing. The three classes
are homogeneous RAYLEIGH with finite positive class/channel scales. All three
relation categories (DIRECT, NEIGHBOR, PRIOR), positive and nonpositive
corrected intensities, six-face border preservation, smoothing, and nonzero
inverse-shape prior are exercised.

Regenerate the NIfTI inputs with:

```sh
python scripts/generate_toads_rayleigh_membership_fixture.py
```

Regenerate the Java composition reference with:

```sh
scripts/generate_java_toads_rayleigh_membership_reference.sh
```

The harness pins the selected `MultipleToadSegmentation.java` source SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884` and
matches energy maps, updated memberships, and objective bits for all four
field/PSF combinations. It covers only this fixed-input homogeneous class
pass; centroid/statistics estimation, topology, preprocessing, convergence,
and complete TOADS/MGDM orchestration remain out of scope.
