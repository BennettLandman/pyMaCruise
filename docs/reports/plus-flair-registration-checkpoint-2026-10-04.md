# Restartable plus FLAIR-to-T1 registration checkpoint

P41b wraps the existing ANTsPy rigid-registration adapter in a persisted checkpoint for the selected MaCRUISEplus upstream FLAIR alignment stage. T1 is the fixed/reference image and FLAIR is moving. The stage accepts compressed NIfTI, normalizes inputs to RAS (with the shared warning for non-RAS data), saves `flair_in_t1.nii.gz` on the T1 grid, persists forward/inverse ANTs transforms, and writes an atomic JSON manifest with input/settings signatures and output hashes.

On `resume=True`, the runner skips registration only when the input/settings signature, manifest status, output files and hashes, and warped-image grid/orientation/finiteness all validate. Changed inputs/settings reject resume. A missing or invalid output cache is recomputed. The manifest explicitly records that ANTsPy is a backend substitution and NiftyReg/Aladin numerical equivalence is not established.

## Synthetic validation

Command:

```sh
./.venv/bin/pytest -q tests/test_plus_workflow.py tests/test_registration.py
```

Result: **8 passed in 16.73s** on the pre-commit working tree. The test uses the committed synthetic registration phantom, checks finite RAS output on the T1 grid, persisted transforms and manifest, and monkeypatches registration to prove that valid resume skips the expensive registration call. It also checks that changed settings are rejected. This is focused synthetic integration evidence, not source equivalence, whole-plus validation, or real-data validation. Rerun after the implementation commit before treating this as a completed plan row.
