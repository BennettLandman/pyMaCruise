# Full pytest evidence: LesionTOADS topology NIfTI geometry

- **Implementation commit:** `f3ea4073b4933d0dc3a6a712f3bce59c9ec9763a`
- **Command:** `.venv/bin/pytest -q`
- **Result:** 594 passed, 2 expected RAS-reorientation warnings, 79.50 seconds.
- **Focused selection:** `tests/test_lesion_atlas_priors.py` — 24 passed in 0.76 seconds.

The selected VandyCRUISE descriptor's optional `Dimensions` row is now checked
against the RAS-oriented NIfTI array shape. Malformed dimension lists and
shape mismatches reject with a descriptor error. The bundled atlas resolves
to `(172, 216, 164)`, remains integer-valued, and is RAS-oriented. A synthetic
LPS image still warns and reorients; a minimal descriptor without `Dimensions`
continues to work. Tests also verify that this loader does not require a shape
map to exist.

This validates the compressed-NIfTI asset contract only. The selected Java
`LesionToadAtlas.loadTemplateImage` reads a headerless `Nx*Ny*Nz` byte payload;
the self-contained repository contains the converted `.nii.gz` asset and not
that payload. The Python loader intentionally uses NiBabel and this report
does not claim byte-for-byte comparison against the unavailable raw file or
the Java loader.
