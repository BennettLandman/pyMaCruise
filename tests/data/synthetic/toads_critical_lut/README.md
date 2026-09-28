# TOADS 18/6 critical-point LUT reference

This synthetic fixture records 68 binary 3×3×3 neighborhoods, their integer
keys, and regular/critical decisions produced by the selected VandyCRUISE
`CriticalPointLUT.keyFromPattern(boolean[][][],...)` and `get(int)` methods.
Four edge cases (empty, one corner, two opposite face neighbors, and full)
are followed by 64 deterministic cases generated with Java `Random(186)`.
The center bit is encoded in each fixture pattern to verify that both Java and
Python ignore it when forming the key.

The class copy in `tests/java_reference/upstream/.../CriticalPointLUT.java` is
byte-for-byte copied from the selected source; SHA-256:
`50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f`.
The compressed table is the bundled selected-source resource; SHA-256:
`cbf1e339d78a0a2c673c939ee2c8a30ee3a35bfd9603d07ed2d48864bb0e88be`.
This validates only key encoding and the source LUT's one-neighborhood
regularity predicate. `MultipleToadSegmentation.isHomeomorphicSegmentation`
also checks changes to the old object and neighboring label pairs/triples;
those multi-label constraints are not represented by this fixture.

Regenerate with a JDK installed:

```sh
scripts/generate_java_toads_critical_lut_reference.sh
```
