# LesionTOADS data-energy-to-membership integration

Two selected-Java synthetic voxel cases start from two-channel intensity
samples, tissue centroids, and channel weights. They compute raw unfielded
class/lesion/black-hole residual energies, add six-face smoothing and
shape/similarity prior terms, then convert, distance-modulate, normalize,
and accumulate local energy for LESION and BLACKHOLE branches.

The neighbor stencils and shape/similarity vectors are fixed synthetic
inputs. These cases integrate scalar algorithm stages; they do not compute
those spatial inputs or run a subject volume. Regenerate with
`scripts/generate_java_lesion_toads_data_energy_membership_reference.sh`.
