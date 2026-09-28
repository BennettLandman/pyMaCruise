# TOADS shape-relation phantom

The committed synthetic fixture contains a 7×6×5×4 float32 shape-prior
volume and a matching four-slot integer closest-label volume. Both are RAS
`.nii.gz` on a 1 mm isotropic grid. Its deterministic recipe is in
`scripts/generate_toads_shape_relation_fixture.py`.

The fixture covers source direct slots 0–1, neighbor slots 2–3, prior
fallback, competing direct/neighbor matches, zero shape, lower and upper
factor clipping, and empty label slots (`-1`). `java_reference.json` contains
the source-formula factors and prior-only zero-shape suppression decisions
using `priorFactor=0.75f`. It records SHA-256 values for
`MultipleToadSegmentation.java` and VandyCRUISE `Numerics.java`; the latter's
`bounded` implementation is `max(a,min(x,b))`.

Regenerate the NIfTIs with `scripts/generate_toads_shape_relation_fixture.py`
and the Java formula reference with
`scripts/generate_java_toads_shape_relation_reference.sh`.
