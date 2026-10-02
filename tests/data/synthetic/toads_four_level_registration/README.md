# Four-level synthetic TOADS registration fixture

This fixture exercises VandyCRUISE's coarse-to-fine rigid atlas registration
with the source schedule at scales 8, 4, and 2. It is entirely synthetic.

`t1_ras.nii.gz` and `brain_mask_ras.nii.gz` are a 32³ ellipsoidal phantom on a
1 mm RAS grid. `atlas_topology_ras.nii.gz` and the four-channel
`atlas_shape_priors_ras.nii.gz` are 16³; this is the smallest even atlas grid
that retains two samples per axis at scale 8. `generate_fixture.py` recreates
all four NIfTI files deterministically.

No four-level Java result is claimed for this registration run. The selected
Java pyramid harness covers source-faithful subsampling and rigid optimizer
updates at scales 4 and 2 on its own committed fixture; the focused test here
checks the additional scale-8 stage and transform handoff. Scale-specific Java comparisons remain in the existing
`toads_registration_pyramid_optimizer` tests.

## Scale-8 behavior

With one bounded optimizer iteration, the scale-8 result has an unchanged
finite transform but `NaN` energy and stop reason `nonfinite_relative_change`,
after 11 rejected line-search trials. Scales 4 and 2 then run from that
unchanged transform and produce finite results. This records the selected
source-shaped state transition; it does not claim that the compact phantom
achieves an atlas fit. The committed 19×17×15 optimizer atlas cannot test
scale 8 because its coarse z dimension is one, and the common 15×13×11 atlas
is smaller still. This 16³ synthetic atlas retains 2³ samples at scale 8 and
makes the behavior reproducible without subject data.
