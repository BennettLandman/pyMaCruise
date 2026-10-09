# Integrated base-workflow MaACE comparison against Java

## Result

The selected Java MaACE wrapper was run on the exact RAS-cropped gray- and
white-matter inputs emitted by the corrected synthetic base workflow's
CRUISE-core preparation, using the workflow's configured **one diffusion
iteration**. All three terminal products match the saved Python workflow
outputs voxel-for-voxel:

| Product | Changed voxels | Dice / values |
|---|---:|---:|
| Enhanced gray matter | 0 / 5,317,312 | Exact; MAE 0 |
| ACE skeleton | 0 / 5,317,312 | Exact; Dice 1.0; 70,225 foreground voxels |
| Thinned skeleton | 0 / 5,317,312 | Exact; Dice 1.0; 47,081 foreground voxels |

To reproduce the workflow's iteration cap without editing the supplied Java
source, the reference generator compiles a **build-only copy** of
`MaAnatomicallyConsistentEnhance.java` with its `max_iter` constant changed
from 25 to 1. The original source remains unchanged; the generated
`spec.json` records the override. The change matches the synthetic workflow
command, which uses `--ace-diffusion-iterations 1` to keep the integration run
small. The prior exploratory 25-iteration comparison on these one-iteration
workflow outputs was a parameter mismatch and is superseded by this result.

The Java-reference diffusion stage was also compared with Python's
source-faithful `ace_anisotropic_diffusion(..., source_faithful=True)` on the
same held-fixed inputs. Both had **0 NaN voxels**. Among 5,317,312 finite
voxels, 31,075 values differ by float roundoff; MAE is `6.34e-11` and maximum
absolute error is `3.65e-7`. This reports the historical zero-tensor behavior
explicitly rather than folding it into a final-output tolerance.

This is a selected-source stage comparison on synthetic data, not full
MaCRUISE composite equivalence or real-MRI validation.

## Reproduction

```sh
BASE=build/base-workflow-2026-10-09-retained-center
./.venv/bin/python scripts/generate_java_ace_reference.py \
  --place-to-add-csf 0 --diffusion-iterations 1 \
  --input-root "$BASE/java-ace/input" \
  --output-dir "$BASE/java-ace/current-source-faithful-reference"
./.venv/bin/python scripts/compare_ace_current_stage.py \
  "$BASE/java-ace/current-source-faithful-reference" \
  "$BASE/stages/cruise-core/attempt-01" \
  --diffusion-iterations 1 \
  --report "$BASE/java-ace/current-comparison.json"
```

The exact output metrics are in
[`ace-integrated-base-source-comparison-2026-10-09.json`](ace-integrated-base-source-comparison-2026-10-09.json).
The corresponding direct Java replay used `-Xmx2g` and remained below the
8 GB watchdog (1.037 GB peak sampled RSS). Inputs, reference products, and
intermediate Java fields are generated under ignored `build/` paths.
