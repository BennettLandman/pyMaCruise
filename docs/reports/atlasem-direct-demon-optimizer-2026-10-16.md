# AtlasEM direct Demon rigid optimizer comparison — 2026-10-16

## Finding

The post-recrop AtlasEM Python registration now matches the supplied
`DemonToadDeformableAtlas` class on a controlled synthetic input. The class
was compiled and executed directly, rather than modeled with the related TOADS
implementation. On both a one-iteration run and a three-iteration cap, all
six final float32 transform parameters match bit for bit.

This exposed a source behavior difference. The Java
`computeRegistrationCoefficients()` loop samples every image voxel selected
by its stride-two scan and does not test class shape bounds. It recomputes
shape bounds after registration. The Python adapter had applied
pre-registration bounds during optimization. The AtlasEM-specific adapter
now supplies full image bounds to the shared kernel; other TOADS registration
paths retain their selected class-bound semantics.

## Synthetic case and reference

- Image and shape grid: `13 × 11 × 9` voxels.
- Image spacing: `(1.0, 1.2, 1.5)` mm; atlas spacing: `(0.8, 1.0, 1.25)` mm.
- Two classes: background and one registered smooth shape; membership values
  form an offset smooth target.
- Initial transform: six zeros; source mode: single rigid, matrix-first,
  stride two, offset zero, source-faithful coefficient order.
- The harness loads the same float32 arrays in the Java class's x/y/z order,
  initializes its public class/image setters, selects the registered class,
  and executes `initShapeRegistration()` followed by `registerShapes()`.

The selected source hashes are pinned in the reference runner. The primary
class hash is `96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a`.
The runner uses a small headless `ArrayUtil` conversion helper because that
JIST utility is not in the supplied checkout; it does not replace optimizer
or interpolation routines. MIPAV display calls use the existing quiet test
stub.

| Iteration cap | Java transform float32 bit patterns | Python transform bit patterns | Result |
|---:|---|---|---|
| 1 | `1054641176,1054641176,1054641175,1054641175,1054641175,1054641176` | same | Exact |
| 3 | `1054641176,1054641176,1054641175,1054641175,1054641175,1054641176` | same | Exact |

## Validation

Focused command:

```sh
PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_demon_rigid_java_parity.py \
  tests/test_atlas_em_registration_from_descriptor.py --tb=short
```

Result: **6 passed in 3.76 s**. Full pytest on code commit
`e11753d350adc9db80a2e8619b9280f7f83bb4fc`: **1,648 passed, 4 skipped,
0 failed, 21 expected orientation warnings in 736.03 s**. Peak sampled RSS was
**14.206 GB** under a 16 GB watchdog; the watchdog did not terminate the run.
The exact full-suite command and result are in
[`full-pytest-2026-10-16-atlasem-direct-demon.md`](full-pytest-2026-10-16-atlasem-direct-demon.md).

This validates the selected class's synthetic single-scale optimizer output.
It does not validate its initial pyramid optimizer, selected `ibsr8obj`
AtlasEM assets, human AtlasEM convergence, pre-LesionTOADS integration, or
whole-PLUS equivalence.
