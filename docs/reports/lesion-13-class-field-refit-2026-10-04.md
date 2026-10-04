# LesionTOADS 13-class exact field refit — 2026-10-04

## Ported behavior

The exact-membership global Chebyshev field fitter and its channel-ordered
adapter previously required the selected 11-class base MaCRUISE descriptor.
They now accept an explicit descriptor-ordered class-type sequence while
retaining the original 11-class layout as the default. The correction-weight
dispatch follows the selected Java method: optimized_csf, optimized_gm, and
optimized_wm classes use their configured power table; regular obj and lesion
classes receive unit weight. The selected LesionTOADS wrapper leaves its
optimized exponents at 1.0.

The committed synthetic fixture
tests/data/synthetic/lesion_inhomogeneity_13class_refit/ adds the selected
13-class profile to an existing small exact-field phantom. It includes a
lesion class and a thirteenth object class, their descriptor type order,
two-channel-ready class centers, and RAS .nii.gz inputs. Its generator and
derivation are documented alongside the fixture.

## Validation

The updated formula harness accepts explicit descriptor class types and
compares the Python 13-class normal matrix and right-hand side exactly against
the selected Java source formula. It also checks the fitted coefficients and
field values within the existing float32 tolerance, for both channels in
source order. The baseline 11-class degree-2/3 and multichannel tests remain
in the same run.

Focused command:

    ./.venv/bin/pytest -q tests/test_toads_exact_polynomial_field_degrees23.py tests/test_toads_multichannel_exact_field_refit.py tests/test_toads_post_thinning_exact_field.py tests/test_toads_post_thinning_exact_field_ramp.py tests/test_lesion_inhomogeneity_transition.py tests/test_toads_inhomogeneity_fit_degree1.py

Result on commit a463d83: **24 passed, 0 failed, 0 skipped in 10.75s**.
This is synthetic selected-formula evidence. No full pytest suite was run for
this slice. The refit kernel is now usable with the 13-class descriptor, but
the full LesionTOADS loop does not yet assemble or consume this refit in its
production workflow; P39 and the plus workflow remain incomplete.

## Source anchors

- InhomogeneityCorrection2.computeExact3DPolynomialField and its optimized
  class weight dispatch, source SHA-256
  6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8.
- AlgorithmLesionToads.runAlgorithm configures channel-wise correction,
  source SHA-256
  556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
- Modified formula harness SHA-256
  d74e2719b47742ecf7eacf76b26d072c67b6c42ab2b50c0fb51d6d6661c85f24.
