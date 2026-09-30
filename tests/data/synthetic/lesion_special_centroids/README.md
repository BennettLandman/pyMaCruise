# LesionTOADS LESION/BLACKHOLE centroid phantom

`special_centroid_inputs.nii.gz` is a deterministic 5×4×3, 1-mm RAS phantom
with 12 float32 frames: labels, five class memberships, two image channels,
two fields, and lesion/black-hole membership maps. The class types are
OBJECT, LESION, BLACKHOLE, OBJECT, OBJECT; the fifth class is absent.
Deterministic map overrides force class/lesion, lesion/black-hole, and
all-three ties, each winner branch, and an exact zero lesion residual.

Generate the NIfTI with `scripts/generate_lesion_special_centroids_fixture.py`
and the selected Java reference with
`scripts/generate_java_lesion_auxiliary_centroids_reference.sh`. The harness
compiles against the preserved upstream `PowerTable.java`; it records
branch/tie counts and raw float32 class, lesion, and black-hole centroid bits
for `none`, `linked`, `damped`, `prior`, and fallback modes, with and without
fields. This validates centroid routing on supplied memberships only; it does
not estimate lesion/black-hole maps or run the full LesionTOADS iteration.
