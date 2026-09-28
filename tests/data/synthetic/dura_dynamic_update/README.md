# Dura dynamic-pressure one-point Java reference

This fixture records two direct calls to the selected VandyCRUISE
`DynamicPressureTGDM.updatePoint` implementation: one interior voxel with
positive pressure and one boundary voxel with negative pressure. Inputs are
deterministic 5×5×5 scalar fields defined by the formulae and parameters in
`java_reference.json`. The reference gradient is calculated by the selected
`FastMarchingGradient.gradient(..., MAGNITUDE)` implementation through the
actual `DynamicPressureTGDM.setReferenceLevelSet` method. The selected
`DynamicPressureCortexProfile.getEffectiveForce` is also used.

The repository copy of `DynamicPressureTGDM.java` is byte-for-byte identical
to the selected source (SHA-256
`53b872a5828539bb51d70db0c048d1fd8162f6b7a6b63176266ef7e2e796b6fc`). The
profile copy is `f85912d621198de297e65d6b9205505b3ac2de72e67e6840f3968c58a850f505`;
the gradient source is `cf857c0c62e2c7a219e9e530c863b4e40ae594b55e90a355cd7fdf6e0709d3a3`.
The selected `GenericTGDM.java` source hash is
`70f5b21541e8a6c6139f1bd227c4ece12e8818db717f2a78ef50d69fe89f08d1`.

## Validation boundary

The dynamic subclass overrides GenericTGDM's point update, so this harness
executes that real override, its reference-gradient setup, and its real dynamic
profile. A test-only `GenericTGDM` state seam supplies the inherited fields and
`NBPoint` type needed to make the override callable. It does **not** compile or
exercise VandyCRUISE `GenericTGDM`'s evolution loop, topology veto, narrow-band
schedule, or distance reinitialization. Compiling the entire base class pulls
in unrelated JIST mesh/spherical-map and MIPAV fast-marching classes; those
dependencies are not part of this focused executable reference. Therefore the
JSON is one-point force/update evidence only, not a dura-removal comparison.

## Regenerate

From the repository root, with a JDK installed:

```sh
scripts/generate_java_dynamic_pressure_update_reference.sh
```

The harness has no imaging input. Its arrays use x/y/z index order and are
only numerical neighborhoods; it does not write NIfTI files or represent
subject data.
