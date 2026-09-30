# LesionTOADS Sulcal-CSF/MPRAGE data-term reference

The Java formula fixture contains ten ordinary-object cases for selected
`Sulcal-CSF` data-energy behavior. In the unfielded source, the special T1
direct squared residual or MPRAGE PowerTable term is followed by a standalone
generic weighted all-channel loop: the apparent `else` before that loop is
inside a comment. In the fielded source, the active `else` makes the exact
Sulcal-CSF branch T1-only. Cases also cover field multiplication, fallback to
generic weighted multichannel energy when T1 is absent, exact-name gating,
and a rounding-sensitive check of source accumulation order.

Regenerate with
`scripts/generate_java_lesion_toads_sulcal_csf_reference.sh`. The fixture
validates only this class data term; spatial smoothing, priors, inversion,
and membership normalization are excluded.
