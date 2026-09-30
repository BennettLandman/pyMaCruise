# P39zzzzv continuing-iteration phantom

The NIfTI file is a deterministic 9×9×9 folded, asymmetric, multi-channel RAS
phantom. Frames 0–1 are T1/T2-like image channels; 2–3 are fixed intensity
fields; 4 is the initial segmentation; 5–8 are class memberships; 9–11 are
lesion, black-hole, and outlier memberships; 12–15 are fixed shape maps; 16
is the starting label state; 17 is the initial skeleton; and 18 is ordering.
The generator is in this directory.

The test starts with the actual Python outputs of the validated P39zzzzs,
P39zzzzt, and P39zzzzu components, then composes the next thinning and
membership update. The relation map created after prior growth is intentionally
retained through thinning, matching the selected wrapper; per-voxel direct vs.
prior membership dispatch is checked independently against the selected Java
lookup rule. Correction fields and atlas maps remain supplied fixed inputs.

Pinned source provenance: `AlgorithmLesionToads.java` SHA-256
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`,
`LesionToadSegmentation.java` SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`, and
`LesionToadAtlas.java` SHA-256
`718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e`. The
Java formula excerpt, reference JSON, hashes, and fixture are committed here;
running tests does not require the original checkout.
