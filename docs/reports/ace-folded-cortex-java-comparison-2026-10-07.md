# Folded-cortex ACE comparison against Java — 2026-10-07

## Scope

This is a selected-source comparison of the MaCRUISE ACE wrapper on the
committed 41×41×41 folded cortex phantom. The comparison uses the exact
committed `gray_matter.nii.gz`, `white_matter.nii.gz`, and
`multi_atlas_csf.nii.gz` arrays as inputs to both implementations. The Java
wrapper is configured with the base MaCRUISE setting, `multiAtlasCSF` merged
after diffusion. This is synthetic evidence for the selected ACE profile, not
a full MaCRUISE run or real-MRI validation.

Java source hashes:

- `MaAnatomicallyConsistentEnhance.java`: `f3cd1c9f4e90cd58ae5cdf37237d90f5a57c59f1a6d6870c9106f4850174cffb`
- `AnatomicallyConsistentEnhancement.java`: `706e20dba6eb58584d0026b34f2cf4965fde532fc88f5304488dd0d8baee5f46`
- `ThinACE.java`: `87c4c99ae01e2e3fb2ab744b71c1dc9367e6c625969badf3ab51158fe5b78203`

The generated NIfTI references and per-file SHA-256 manifest are in
[`tests/data/synthetic/folded_cortex/ace_java_reference`](../../tests/data/synthetic/folded_cortex/ace_java_reference).

## Reproduction

```sh
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 -- \
  .venv/bin/python scripts/generate_java_ace_reference.py \
  --input-root tests/data/synthetic/folded_cortex \
  --output-dir tests/data/synthetic/folded_cortex/ace_java_reference

.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 -- \
  .venv/bin/pytest -q tests/test_ace_reference.py
```

The Java process is configured with a 2 GiB maximum heap. The focused Python
ACE reference suite passed **3 tests in 10.04 s** on committed implementation
`f00117f`, with sampled peak RSS
**0.156 GB** under the 4 GB watchdog. This is focused evidence only; no full
repository suite was run for this slice.

## Findings

- Seed, morphology, and diffusion membership stages that are deterministic on
  this input match their selected Java references exactly. Other early float
  stages match within `5.73×10⁻⁶` maximum absolute error and `7.22×10⁻⁷` mean
  absolute error.
- The narrow-band evolved level set differs by at most `0.011816` intensity
  units (mean `1.244×10⁻⁵`). Its thresholded mask Dice is **1.0**; symmetric
  mean surface distance is **1.15×10⁻⁸ mm** and Hausdorff is **1.88×10⁻⁶ mm**.
- Fast-marching weighted distance differs by at most `0.020771` (mean
  `4.03×10⁻⁶`). The pre-cleaned and cleaned scalar CSF membership arrays differ
  at **3 of 68,921 voxels**, by one intensity unit; their nonzero supports are
  identical (Dice **1.0**).
- The selected wrapper's enhanced-GM, CSF-skeleton, and thinned-skeleton NIfTI
  outputs match Java voxel-for-voxel. Both skeleton masks have Dice **1.0** and
  identical measured surfaces (zero symmetric mean and Hausdorff distance).

These results meet the project's segmentation and surface tolerances on this
folded fixture. They also preserve the small scalar-membership deviation in
the evidence instead of treating the intermediate floating-point outputs as
bit-identical. The 25³ concentric reference remains covered by the same focused
suite.
