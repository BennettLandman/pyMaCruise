# Default TOADS atlas dispatch fixture

This committed synthetic fixture is a folded, asymmetric 9×8×7, two-channel
phantom with true-RAS geometry and anisotropic millimeter spacing. Its class
axis follows the selected default atlas descriptor exactly: MASK background,
two OPTIMIZED_CSF, five OPTIMIZED_GM, then three OPTIMIZED_WM. Distinct
per-class centroids and variance weights exercise repeated category dispatch;
all instances of a category share its channel-specific PowerTable exponents.
The fixture includes active shape maps, four closest-label slots, field and
PSF-bias maps, and initial memberships.

Regenerate the inputs with:

```sh
python scripts/generate_toads_default_atlas_fixture.py
```

Regenerate the selected Java reference with:

```sh
scripts/generate_java_toads_default_atlas_reference.sh
```

The Java harness pins `MultipleToadSegmentation.java` SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884` and
`PowerTable.java` SHA-256
`5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678`.
Reference cases cover field/PSF and smoothing/prior on/off, including the
source behavior that prior-only classes are suppressed when the prior
coefficient is zero. The harness validates one fixed-input exact-membership
dispatch slice, not statistics estimation or the complete TOADS workflow.
