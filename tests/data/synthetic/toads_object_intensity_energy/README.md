# TOADS OBJECT intensity-energy fixture

The fixture records the ordinary `OBJECT` branch of selected
`computeExactMemberships` for a synthetic 4×3×2 grid, three modalities, and
four classes. It covers uncorrected and field-corrected intensities, with and
without an explicit PSF-bias array. Inputs and outputs are stored as float32
bits so the Java/Python comparison is exact. There is no medical image data.

Regenerate with
`scripts/generate_java_toads_object_intensity_energy_reference.sh`. The
selected source file is `MultipleToadSegmentation.java`, SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
