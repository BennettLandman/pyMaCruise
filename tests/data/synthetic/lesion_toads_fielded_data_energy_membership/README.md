# LesionTOADS field-corrected data-energy-to-membership integration

Two selected-Java synthetic voxel cases start from two-channel intensity and
nonunit field samples, tissue centroids, and channel weights. They multiply
each channel by its field sample before computing class/lesion/black-hole
residual energies, then add six-face smoothing and shape/similarity prior
terms, convert, distance-modulate, normalize, and accumulate local energy
for LESION and BLACKHOLE branches. The Python result is compared bitwise with
the Java float32 harness.

The neighbor stencils and shape/similarity vectors are fixed synthetic
inputs. These cases integrate scalar algorithm stages; they do not compute
those spatial inputs or run a subject volume. Regenerate with
`scripts/generate_java_lesion_toads_fielded_data_energy_membership_reference.sh`.
