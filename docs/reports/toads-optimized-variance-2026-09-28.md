# P18am: optimized-tissue image variance weights

## Scope

This source slice implements `MultipleToadSegmentation.computeExactImageVariance`
for optimized CSF, GM, and WM classes. For every modality, class, and voxel in
the full image grid (including borders), it accumulates the membership-squared
weighted fourth-power residual and the corresponding stabilized membership
weight. It then produces normalized inverse-variance modality weights and
the square root of each unnormalized inverse variance. An optional
multiplicative field is applied to image intensities before residual
calculation.

The implementation assumes every supplied membership class is an optimized
tissue class. It does not yet cover OBJECT, OUTLIER, or LESION branches, and
does not claim a complete TOADS update or segmentation.

## Synthetic reference

The 7×6×5×3 1 mm RAS image and field volumes are reused from
`tests/data/synthetic/toads_optimized_likelihood/`; the first three
membership channels are reused from
`tests/data/synthetic/toads_direct_membership/`. Three centroid rows from the
selected optimized-tissue initialization are recorded in the Java harness.
No duplicate image fixtures were added.

The standalone Java harness transcribes the optimized-class branch of
`computeExactImageVariance` and records selected source SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884` in
`tests/data/synthetic/toads_optimized_variance/java_reference.json`.
Regenerate the JSON with
`scripts/generate_java_toads_optimized_variance_reference.sh`.

## Validation

`tests/test_toads_optimized_variance.py` compares all three normalized
inverse-variance values and image-standard-deviation values exactly as
float32, both with and without the field. It also verifies that a single
modality normalizes to 1 and that zero residuals produce the source's
`INF=1e10f` per modality before normalization (0.5 each for two modalities).
The focused module passes **3 tests**.

The full configured suite passes **326 tests in 51.59 seconds** on
implementation commit
`6862499021ca1b348559b0cf6f4730a7e4fe6a3a`; the exact run is recorded in
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md). This establishes
formula agreement for the tested synthetic optimized-tissue slice only.
