# MaCRUISE ACE synthetic Java reference

The compressed NIfTI outputs here were produced by the supplied MaCRUISE
`MaAnatomicallyConsistentEnhance.java` and VandyCRUISE ACE/volume Java sources
using the committed 24³ concentric tissue phantom in the parent directory.
The source runner uses the selected V3 configuration: threshold 126.9,
Gaussian alpha 1.8, tensor alpha 1.6, 25 diffusion iterations, updated
multi-atlas CSF merged after diffusion, and the source narrow-band settings.

Outputs are enhanced gray matter, the cleaned ACE skeleton before thinning,
and the VandyCRUISE ThinACE result. All images have a 1 mm isotropic RAS
affine. Source inventory, provenance, and component comparisons are recorded
in `tests/java_reference/README.md` and `docs/ace-port.md`.

Regenerate using `python scripts/generate_java_ace_reference.py`. A JDK and
the project's Python dependencies are required. Compiled classes and binary
interchange files stay in ignored `build/`.
