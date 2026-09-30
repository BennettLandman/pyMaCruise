# TOADS MASK intensity-energy fixture

Selected Java output for the `MASK` intensity branch and its
`maskDistance(image, 0)` helper on a synthetic 3×2×4 volume with three
modalities and four classes. Inputs include positive, negative, exact-zero,
and signed-zero samples with nonuniform float32 weights. Every source and
output value is recorded as float bits. There is no medical image data.

Regenerate with
`scripts/generate_java_toads_mask_intensity_energy_reference.sh`. The source
is `MultipleToadSegmentation.java`, SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
