# LesionTOADS voxel-membership normalization reference

This synthetic selected-Java formula fixture contains eight one-voxel cases
for the tail of `computeMemberships`, starting from class and auxiliary
memberships and class/auxiliary energies already computed by earlier steps.
It records source-order denominator accumulation, normalized class,
lesion, black-hole, and outlier memberships, and the resulting local energy.

Cases cover ordinary classes, LESION, BLACKHOLE, OUTLIER/OUTMASK, a
nonpositive denominator, repeated LESION auxiliary division, optimized
lesion, and float32-sensitive sums. In particular, the harness retains the
source behavior that a shared auxiliary map is divided once per matching
class type. Regenerate with
`scripts/generate_java_lesion_toads_voxel_membership_normalization_reference.sh`.
This is a scalar formula harness, not full voxelwise membership estimation.
