# LesionTOADS Sulcal-CSF/MPRAGE data-term reference

The Java formula fixture contains nine ordinary-object cases for selected
`Sulcal-CSF` data-energy behavior. It covers the T1-only direct squared
residual, T1-only MPRAGE PowerTable lookup, field multiplication, fallback to
generic weighted multichannel energy when T1 is absent, and exact-name gating
for the Sulcal-CSF special case.

Regenerate with
`scripts/generate_java_lesion_toads_sulcal_csf_reference.sh`. The fixture
validates only this class data term; spatial smoothing, priors, inversion,
and membership normalization are excluded.
