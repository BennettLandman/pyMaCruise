# LesionTOADS capped distance support phantom

This fixture validates only the small distance-to-mask support method from
the selected VandyCRUISE source, `LesionToadSegmentation.fastDistance` with a
boolean mask and byte output. It does **not** validate LesionTOADS
segmentation equivalence, class priors, atlas registration, or lesion labels.

The synthetic cases cover an empty mask, a full mask, an object at the image
boundary, an isolated object with caps 1, 5, and 9, and diagonal-only offsets.
Arrays follow Java's `[x][y][z]` indexing; flattened reference values traverse
x outermost, then y, then z. Expected results are integers and are compared
voxel-by-voxel, including the source sentinel `-1`.

The source algorithm reports 0 inside the object, six-connected voxel-step
distances outside, and `-1` where the cap is exceeded or no object exists. It
ignores physical voxel spacing. Its post-scan cap check has a small edge case:
with `max_distance=1`, the mandatory first propagation pass also assigns
distance 2. The Python port preserves this result. It accepts positive caps
from 1 through 127, the useful positive range of the source's signed Java byte;
invalid zero, negative, or overflowing caps are rejected rather than
reproducing Java byte wraparound behavior.

Regenerate the Java reference JSON with:

```sh
javac -d /tmp/pymacruise-p39-reference tests/java_reference/LesionToadsDistanceReferenceHarness.java
java -cp /tmp/pymacruise-p39-reference LesionToadsDistanceReferenceHarness \
  > tests/data/synthetic/lesion_toads_distance/java_reference.json
```

The test harness preserves the selected method body semantics in a standalone
Java class so the primitive can be checked without MIPAV/JIST runtime setup.
This is a support-kernel comparison, not a full invocation of the original
`LesionToadSegmentation` class.
