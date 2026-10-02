# LesionTOADS output module CLI

The selected terminal-output arrays could already be written through the
Python API, but there was no command-line adapter for a caller that has those
arrays as NIfTI files. Added `python -m pymacruise.lesion_output_cli` to load
the output arrays and original T1 reference, route the selected output mode,
and delegate RAS orientation, type conversion, source-derived filenames,
collision checks, and `.nii.gz` serialization to the existing P39 output
writer.

This is an output-materialization command only. It does not run LesionTOADS,
infer segmentation results, resample input product arrays, or claim workflow
equivalence. The caller must provide already-computed arrays in RAS voxel
order on the original T1 grid. `--mode-data` is `[x,y,z,channels]` for
`class+member`, and `[4,x,y,z]` in CSF/GM/WM-fill/WM-mask order for the CRUISE
and Dura branches. The reference is normalized to RAS by the shared output
writer, which emits its standard warning when reorientation is needed.

Synthetic CLI validation used an asymmetric 5×4×3 LAS reference and
four-channel membership fixture. The subprocess produced the three expected
compressed outputs, preserved membership voxel values/channel order and
float32 type, and wrote RAS-oriented images. Focused validation:

```text
./.venv/bin/pytest -q tests/test_lesion_output_cli.py tests/test_lesion_output_io.py
11 passed in 1.52s
```

Invocation help is available with:

```text
python -m pymacruise.lesion_output_cli --help
pymacruise-lesion-output --help
```

`pyproject.toml` now registers `pymacruise-lesion-output` as an installed
console script. The editable local installation was refreshed and the
installed entry point's `--help` output was verified. A fresh isolated wheel
build/install smoke test has not yet been run, so wheel packaging remains
unverified.
