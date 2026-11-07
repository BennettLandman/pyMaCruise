# BASE axial windowed-sinc reference

Synthetic 13×14×11 masked intensity field with an impulse, a low outlier and
anisotropic float32 spacings. `java_expected_ras.npz` contains the independently
computed selected Java `ReorientVolume` AXIAL/UNCHANGED/WINDOWED_SINC output,
converted from the SPECTRE LPS buffer to RAS by flips only. The Java transform
uses MIPAV `AlgorithmTransform`/`AlgorithmWSinc` with clipping disabled.

Generated with `scripts/compare_human_spectre_reorient.py`; the source SHA-256,
operation and fixture identities are in `provenance.json`. No patient data.
The Python test reads the saved Java values and does not regenerate its oracle.
