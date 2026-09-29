# LesionTOADS lesion/black-hole and OUTMASK integration

Two selected-Java synthetic voxel cases combine a LESION or BLACKHOLE
auxiliary map with an OUTMASK/outlier mixture, then pass all memberships
through the shared denominator, normalization, and local energy expression.
Cases cover GM-zero plus brainstem for LESION and ventricle precedence for
BLACKHOLE (where the source leaves its black-hole map unchanged), alongside
prior/direct OUTMASK mixing.

Inputs are already converted class/auxiliary memberships and energies. This
fixture validates cross-auxiliary denominator interactions, not energy-term
construction or volume-wide segmentation. Regenerate with
`scripts/generate_java_lesion_toads_auxiliary_membership_integration_reference.sh`.
