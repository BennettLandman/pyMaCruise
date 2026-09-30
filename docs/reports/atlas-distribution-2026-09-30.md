# Atlas distribution verification — 2026-09-30

This is a packaging check for plan row P42b, not a full test-suite report and
not evidence that P42b is complete.

## Run

- Command: `.venv/bin/python scripts/verify_atlas_distribution.py`
- Script result: success (exit status 0).
- Repository `HEAD` at invocation: `06e3c821ee731c63487179e6095f317199028415`.
- The working tree contained uncommitted source and documentation changes. The
  script built a source copy from the checkout at invocation; therefore these
  hashes are artifact evidence for that tree, not a clean-commit build.
- Build artifacts were preserved beneath
  `tmp/build/atlas-distribution-45403-1j8wj9zx/`.

## Results

- All 77 repository atlas files appeared exactly once in both the source
  distribution and wheel, with hashes matching their source files.
- The isolated wheel install resolved the bundled default TOADS (11 classes),
  TOADS 2012 (13 classes), and LesionTOADS (13 classes) inputs without the
  checkout on `PYTHONPATH`.
- All 35 packaged NIfTI files were `.nii.gz` and reported RAS orientation.
- Representative atlas references and the 172×216×164 topology image loaded
  from the installed wheel's package resources.
- Source distribution: 140,340,642 bytes;
  SHA-256 `f61cff6a04aa373b307310d9a9ac2455b4be45b8ff650c708b06248d8e2189e7`.
- Wheel: 140,336,367 bytes;
  SHA-256 `11db1f444895f4b4d11aa8b083d3024fa297c24cd418b0bb9799328eb30eb9b6`.
