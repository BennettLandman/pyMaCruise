# LesionTOADS ordinary classification update phantom

This synthetic fixture pins the selected VandyCRUISE
`LesionToadSegmentation.updateClassificationFromMemberships()` method
(source SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`;
wrapper SHA-256 `556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`).
The Java harness is a selected-method excerpt, not a full LesionTOADS run.

The `multiple_classes` fixture covers first-wins ordinary ties, all-zero
memberships, lesion strictly winning and tying the ordinary winner, ordinary
class-zero behavior, template byte labels that differ from membership indices,
signed byte labels, prior classification values, and the returned changed-voxel
count. The `one_class` fixture records the source's special result when the
ordinary class axis has only class zero. Memberships and lesion memberships are
finite float32 values; classification and template labels are signed bytes.

All stored images are `.nii.gz` NIfTI with a one-millimeter identity RAS affine.
Memberships store the class axis last. Regenerate the pinned Java JSON with
`sh scripts/generate_java_lesion_classification_from_memberships_reference.sh`
and the NIfTI fixtures with
`python tests/data/synthetic/lesion_classification_membership_update/generate_fixtures.py`.

This small routine updates an ordinary classification byte array in place. It
does not implement or imply the full LesionTOADS workflow.
