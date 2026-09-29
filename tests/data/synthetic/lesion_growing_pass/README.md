# LesionTOADS growing queue exhaustion phantom

`folded_labels.nii.gz` is a deterministic 9×9×9 int16 RAS label image with
1.5×2.0×2.5 mm voxel spacing. It contains a folded stair-stepped class-1
band, a compact class-2 region, and a narrow class-3 band. There are no
subject images. The two selected Java reference runs use this same label
geometry with deterministic formula-generated float32 memberships, ordering
values, class volumes, and three skeleton seeds.

The additive and max lesion-score runs exhaust the min-heap and record each
pop, topology outcome, neighbor insertion, and final segmentation, skeleton,
availability, ordering bits, volume bits, and queue state. This validates
the selected growing pass using supplied memberships. It does not cover
membership estimation, competition, thinning, or complete LesionTOADS output.

Regenerate the label image with `make_phantom.py` and the Java references with
`scripts/generate_java_lesion_growing_pass_reference.sh`.
