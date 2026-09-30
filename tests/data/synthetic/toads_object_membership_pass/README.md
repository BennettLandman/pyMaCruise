# TOADS ordinary OBJECT membership-pass fixture

This fixture composes the selected exact-membership ordinary-OBJECT path on a
synthetic 4×4×4 volume with two modalities and three classes. It records a raw
case and a field-corrected case with PSF bias, along with initial memberships,
transformed shape values, closest-label support, class/model weights,
smoothing/prior parameters, final memberships, and objective values. The
ordered sweep updates in place and the output is compared bit-for-bit to the
selected Java formula harness. No medical image data is present.

Regenerate with
`scripts/generate_java_toads_object_membership_pass_reference.sh`. The source
is `MultipleToadSegmentation.java`, SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
