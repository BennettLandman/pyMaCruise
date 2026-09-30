# P18cl, P18cm, and P42 validation — 2026-09-30

## Full-suite evidence

- Command: `.venv/bin/python -m pytest -q`
- Tested implementation commit: `82a0d8f495a837422dc17c2b4a5b30b81f6b325b`
- Result: **640 passed, 2 expected RAS reorientation warnings, 160.07 seconds**.
- The warnings are from `test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`, which deliberately supplies LPS-oriented topology volume and paint mask images and verifies visible conversion to RAS.

Focused evidence on the same tree:

- P18cl: `.venv/bin/python -m pytest -q tests/test_toads_preprocessing.py tests/test_toads_registration_sampling.py` — **30 passed**.
- P18cm: `.venv/bin/python -m pytest -q tests/test_toads_crop.py` — **4 passed**.
- P42: `.venv/bin/python -m pytest -q tests/test_cli_standalone.py tests/test_cli.py` — **4 passed**; all seven standalone commands execute in fresh Python subprocesses.

## P18cl — TOADS orientation initialization under RAS normalization

The selected `MultipleToadPreprocess` axial setup maps the RAS axis codes L2R/P2A/I2S to signed spacings `(-sx, -sy, +sz)` and starts with a zero rigid transform. The selected `ParametricTransform.imageToTemplate` uses center-relative voxel coordinates and those signed spacings. Its callers provide dimensions, spacing, axis codes, and image orientation; they do not provide or use a NIfTI world-origin translation.

`toads_ras_geometry` now normalizes image arrays to RAS, warns before reorienting non-RAS inputs, derives the source signs after reorientation, retains the RAS affine metadata, and rejects oblique/sheared grids the selected axis-aligned source mapping cannot represent. Asymmetric images with LPS, axis-permuted, and oblique affines verify data order, warnings, and rejection behavior. Changing only the NIfTI world origin leaves the source signed-spacing/initial-transform mapping inputs unchanged.

A selected-Java `ParametricTransform` landmark harness and Python's field sampler match three RAS-mapping landmarks exactly as float32 coordinates using signed spacing. The compared selected-source files are `MultipleToadPreprocess.java` SHA-256 `31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7` and `ParametricTransform.java` SHA-256 `22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac`; the pinned Java mapping harness is `tests/java_reference/ToadsOrientationMappingHarness.java`. This is landmark/mapping evidence, not full TOADS atlas registration or segmentation equivalence.

## P18cm — TOADS crop-bound discovery

The selected `MultipleToadPreprocess.findCroppingBoundaries()` body (source lines 302–345) and `computeImageCoordinates()` body (lines 653–660) are copied into a small Java method-excerpt harness. The full upstream class cannot be compiled from the supplied checkout because MIPAV/JIST application dependencies are absent; the fixture README identifies this boundary and lists selected VandyCRUISE source hashes. SHA-256 pins cover the excerpt harness and generated reference.

On the committed asymmetric 7×8×6 synthetic RAS case—with anisotropic spacing, canonical `(-1,-1,+1)` signs, nonzero rotation/translation, thresholded image support, and two non-background topology labels—the Python path matches all six inclusive crop bounds exactly. Existing reference comparisons also cover transformed topology labels, crop/uncrop voxel arrays, crop dimensions, and translation updates. This is one selected-method phantom; it does not execute the entire MIPAV/JIST class or establish general/full TOADS equivalence.

## P42 — standalone command-line interfaces

Documented standalone subcommands are `ace`, `gvf`, `tgdm`, `dura`, `surface-extract`, `thickness-levelset`, and `thickness-surface`. Each command invokes its existing Python API, uses compressed NIfTI for image paths, and writes optional `--provenance RUN.json` sidecars with input/output paths, effective parameters, software version, and RAS spatial convention. Surface data are tagged RAS millimetres in legacy VTK. Full TOADS and LesionTOADS commands are intentionally absent because their end-to-end implementations are not ready; the documentation says so explicitly.

Subprocess tests run all seven commands against committed synthetic phantoms and compare written image arrays, mesh vertices/faces, surface measurements, and summary output against direct public API calls. The provenance test checks exact resolved input/output paths and effective nondefault solver parameters. The CLI tests validate standalone interface behavior; they do not validate complete TOADS, LesionTOADS, or MaCRUISE workflow output.
