# Registered LesionTOADS initial membership phantom

This deterministic 9×9×9, two-channel fixture exercises the selected
LesionTOADS wrapper after registered recropping. `images_ras.nii.gz` contains
synthetic T1_MPRAGE/FLAIR ramps, `constructor_labels_ras.nii.gz` contains
interleaved standard-13 class indices, and `transformed_shapes_ras.nii.gz`
contains positive class shape values. All use an identity RAS-mm affine. The
generator and SHA-256 digests are in this directory.

The initial class, lesion, and black-hole centroids are explicit synthetic
seam values. The composition then creates Java-default zero constructor
memberships, outlier membership 0.5, and lesion relation 1. It follows the
selected wrapper order through relations, statistics, channel weights,
ventricle/GM/brainstem distances, the GM-enabled initial membership pass,
loop-entry state, and the first thinning handoff. A high thinning priority
limit is used to ensure the fixture traverses the queue.

Existing component tests carry selected Java source-excerpt/harness
comparisons. This small fixture validates their Python composition; it is not
a composite Java segmentation reference, a registration fixture, or a
whole-workflow LesionTOADS result. Java zero-deviation handling is explicitly
enabled for the membership call. No clinical image data are included.
