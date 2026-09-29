# LesionTOADS FLAIR lesion-weight PowerTable reference

The committed Java formula fixture exercises six `LESION` data-term cases
with `useLesionWeight` enabled. For each case, the class data energy uses
all weighted channel residuals while the auxiliary lesion energy uses only
the FLAIR-channel squared residual passed to the selected Java `PowerTable`.
Cases include cached-grid inputs, out-of-range direct powers, and both raw and
bias-field-corrected intensities.

Regenerate the fixture with
`scripts/generate_java_lesion_toads_flair_weight_reference.sh`. This validates
the specialized lesion data term only; no smoothing, atlas prior, inversion,
or membership update is performed.
