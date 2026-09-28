# Registration and physical-space transfer

**Status:** P33 synthetic validation complete (2026-09-27). The generic
ANTsPy registration/resampling boundary is implemented. MNI-specific
integration against `average305_t1_tal_lin.nii.gz` and its mask remains
unverified because those reference images are not in `data/atlases/`.

## Legacy operation contract

The selected MATLAB wrappers use NiftyReg Aladin with the template/reference
image as `-ref`, the subject image as `-flo`, rigid-only registration for the
MNI path, and the reference-space MNI mask where supplied. Segmentation
resampling explicitly uses nearest-neighbor (`-inter 0`). Continuous images
and level sets omit the interpolation option; NiftyReg's `reg_resample`
default is cubic spline (`-inter 3`). See the
[NiftyReg `reg_resample` source](https://github.com/SuperElastix/niftyreg/blob/master/reg-apps/reg_resample.cpp).
The MNI wrapper applies the estimated affine to segmentations, level sets,
and scientific meshes.

## Python boundary

`pymacruise.registration` uses ANTsPy for rigid/affine estimation and image
resampling. `register_to_reference(fixed, moving, ...)` maps the legacy
reference/floating pair to ANTsPy fixed/moving arguments, optionally applies
the reference-space mask to all registration stages, and requires an output
prefix so transform files persist at a caller-selected location. The
optimizer, metric, sampling schedule, pyramid, deterministic seed, and
precision are configurable and recorded in `RegistrationConfig`; defaults
are explicit rather than claimed to match Aladin's optimizer numerically.

`register_flair_to_t1(t1, flair, ...)` is the MaCRUISEplus adapter: it pins
the transform to rigid-only, uses T1 as the fixed/reference image, and uses
FLAIR as the moving image. Metric and optimizer choices remain configurable.

`apply_image_transform` writes the moving image on the fixed grid. Continuous
images and level sets use cubic B-spline interpolation as the closest
available ANTsPy counterpart to the NiftyReg default; labels use nearest
neighbor and retain their integer dtype. It rejects non-integer label input
and mismatched geometry where a mask is supplied. The output is RAS-oriented
and preserves the fixed image affine/header geometry.

ANTs/ITK uses LPS physical coordinates internally. PyMaCRUISE converts
NIfTI geometry and RAS point/surface inputs at the API boundary. ANTs image
resampling transforms map fixed-grid sample positions back into the moving
image; applying that same transform list to moving-space points requires the
opposite direction. The point helper therefore inverts affine transforms and
evaluates the composite in ITK stack order. This project path currently
supports affine `.mat` point/surface transforms; nonlinear mesh warping is
outside the selected rigid/affine source contract.

## Synthetic validation

The committed 17³ RAS fixture in
`tests/data/synthetic/registration/` contains matching Gaussian volumes and
label spheres separated by a known 2 mm translation, a fixed-space mask, and
valid binary ANTs affine/inverse transforms generated with ANTsPy. Tests
check:

- known forward and inverse image-resampling directions;
- B-spline intensity and exact nearest-neighbor label transfer;
- RAS point and surface mapping, including topology retention;
- fixed-mask geometry checks and a recovered rigid translation;
- non-RAS warning/reorientation and missing-transform/invalid-label errors.

The same fixture includes a square-root intensity-remapped FLAIR phantom;
Mattes mutual-information rigid registration recovers the known 2 mm shift to
within 1 mm, and the plus adapter rejects affine configuration.

The fixture tests establish the Python API conventions and internally
consistent coordinate transfer. They do not establish Aladin-vs-ANTs
equivalence; that requires approved reference outputs or an available MNI
template/mask plus source-generated transform comparisons. In particular,
the spline implementations and registration metrics/optimizers differ.

Regenerate fixture inputs/transforms with
`python tests/data/synthetic/registration/make_phantom.py`, then verify
`SHA256SUMS` from that fixture directory.
