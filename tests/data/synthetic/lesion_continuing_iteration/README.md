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
