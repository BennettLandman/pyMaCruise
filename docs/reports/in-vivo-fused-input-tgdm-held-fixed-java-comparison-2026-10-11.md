# In-vivo fused-input TGDM held-fixed Java comparison — 2026-10-11

## Result

The saved Python base-core run was compared with the selected MaCRUISE-specific Java `MaGenericTGDM` on the same six fused-input crop maps. Python was not rerun. Java completed under a **4 GiB sampled-RSS limit**, with peak sampled RSS **1.231 GiB**.

The three saved terminal surfaces pass the agreed held-fixed comparison limits:

| Surface | Negative-region Dice | Different mask voxels | Symmetric mean surface distance | Hausdorff | Result |
|---|---:|---:|---:|---:|---|
| Central | 1.000000 | 0 | 0.000000005 mm | 0.000014 mm | Pass |
| Inner | 1.000000 | 0 | 0.000000005 mm | 0.000014 mm | Pass |
| Outer | 0.999973 | 53 | 0.000492 mm | 0.937339 mm | Pass |

Acceptance limits are Dice ≥0.97, symmetric mean surface distance ≤0.01 mm, and Hausdorff ≤1 mm. Surface distances were computed between Lewiner marching-cubes meshes of the held-fixed Java and Python zero level sets, using continuous triangle-surface distances in the shared RAS-mm grid. The Java harness writes an identity-origin NIfTI header; comparison therefore uses the common crop affine from the Python core inputs for both meshes.

Central and inner voxel masks and mesh cardinalities match exactly. The outer mask has 53 differing voxels and Dice 0.9999728. Scalar level-set fields are not bit-identical: their mean absolute differences are 0.02716, 0.01904, and 0.04298 for central, inner, and outer, respectively; maximum absolute differences are 12.94, 12.37, and 14.63. No scalar-field tolerance was specified, so the acceptance result is based on the agreed region and surface metrics.

## Scope and limitations

This is real-human-data, held-fixed source evidence for the TGDM terminal surfaces, not a complete workflow equivalence or clinical validation result. Inputs are the exact fused maps and saved ACE/GVF/core products from the corrected T1-grid human run. The crop is `(141, 173, 131)` at 1 mm isotropic spacing; its inclusive bounds in the full Python core grid are `x=30..170`, `y=58..230`, `z=105..235`.

The harness also writes Java `signed_seed`, `smoothed_initial`, `initial_inner`, `tissue_sign`, and `central_shrink` fields. The saved Python core output does not include corresponding intermediate fields, so this comparison does not claim intermediate-field parity and did not regenerate those Python fields. It compares only `central`, `inner`, and `outer`, the terminal level sets present in both saved runs.

This result validates the selected TGDM calls on the same held-fixed human inputs. It does not assess the upstream TOADS/fusion stages, final MaSeg-refined segmentation, surface extraction equivalence, thickness, or the full base workflow. The surface-distance computation uses marching-cubes meshes from both level-set arrays; it is a TGDM level-set comparison and not a validation of the production JIST surface extractor.

## Selected source and settings

The saved Java harness compiles/calls the MaCRUISE-specific `MaGenericTGDM` from `toads-cruise-scm-2026-09-26/toads-cruise`. Harness source SHA-256 is `ae0606c89e414bb45783be1d08f65e3dd84ae59b40022f54712f77dfe3670ce0`; selected Java engine source SHA-256 is `66124d4180864edc77f6fff280ab817eb0fffc061fdf427c9d8d725f850b211d`.

Settings match the selected wrapper: WM/GM isovalues `0.6/0.4`; initial, first inner, central, central shrink, final inner, and outer counts `6/4/7/1/4/7`; initial level `126.9`; initial maximum distance `8 mm`; no brain-mask field. The MaCRUISE wrapper maps the central iteration count to the outer call.

Input and output hashes, per-field metrics, full affine/shape details, source hashes, and script paths are in the machine-readable report at:

`build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/java-maseg-tgdm-replay/comparison.json`

Human image inputs, Java outputs, and diagnostic scripts remain under ignored `build/`; no image data were committed.

## Reproduction

From the PyMaCRUISE repository root, the held-fixed replay is:

```sh
PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 4 --interval-seconds 0.2 --report-interval-seconds 30 -- \
  java -Xms512m -Xmx3g \
  -cp build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/tgdm-java-source-profile/classes \
  HumanSourceProfileTGDMHarness \
  build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/java-maseg-tgdm-replay/input \
  build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/java-maseg-tgdm-replay/output
```

The saved comparison was computed by loading the Java outputs and completed Python core NIfTIs, cropping the latter with the recorded inclusive bounds, and measuring binary Dice plus continuous triangle-surface distances. The comparator did not invoke the PyMaCRUISE pipeline.
