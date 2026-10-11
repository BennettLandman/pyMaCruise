# P41z next-iteration statistics Java comparison

This case directory documents an input contract; the test creates its small
synthetic NIfTI files in a temporary directory. It does not read ignored
`build/` checkpoints or require the source checkout at test runtime.

## Input contract

`case.json` contains the two channel names, the fixed 13-class object-type
profile, class and auxiliary centroids before refresh, unit tissue power
exponents, lesion power exponents `(0.4, 1.2)`, `norm_type: Identity`, and
`outlier_distance: 0.1`. The six NIfTI files are:

| File | Required data |
| --- | --- |
| `images.nii.gz` | Float32 T1_MPRAGE and FLAIR values, shape `(x,y,z,2)` |
| `memberships.nii.gz` | Float32 class memberships, shape `(x,y,z,13)` |
| `segmentation.nii.gz` | Integer class labels, shape `(x,y,z)`, values 0–12 |
| `lesion_membership.nii.gz` | Float32 lesion membership, shape `(x,y,z)` |
| `black_hole_membership.nii.gz` | Float32 black-hole membership, shape `(x,y,z)` |
| `outlier_membership.nii.gz` | Float32 outlier membership, shape `(x,y,z)` |

All images must already be RAS. The comparison currently fixes the source
profile to `MASK, OPTIMIZED_CSF, OPTIMIZED_CSF, OBJECT × 7, LESION, OBJECT,
OBJECT`, with no correction fields, centroid mode `none`, FLAIR at channel 1,
and Identity channel weighting. The Java Identity branch returns uniform
`1/nc` weights. The three updated statistics stages are compared as raw
float32 bits. NaNs are counted explicitly and never made equal implicitly by
the comparator.

## Scope

The harness is a selected-method excerpt of the VandyCRUISE Java source, not a
full `LesionToadSegmentation` object invocation. The temporary phantom
validates input/output wiring and this selected profile only. A P41z
iteration-6 comparison still needs the exact post-thinning segmentation and
membership state, plus the correct current centroids from its checkpoint.
