# Clean dependency install and synthetic BASE run — 2026-10-15

## Result

The initial attempt in the restricted shell could not resolve the package
index. Retrying with network access succeeded and exposed a separate
compatibility issue: the unconstrained resolver selected SciPy 1.15.3 and
NumPy 2.3.5, but SciPy's `_propack` extension failed to load on this macOS 27
host. SciPy 1.14.1 and NumPy 2.2.6 imported and ran correctly with ANTsPy
0.6.3.

The `ants` extra now constrains NumPy to `<2.3` and SciPy to `<1.15` alongside
ANTsPy 0.6.x. A fresh target installation from the updated project metadata
resolved these versions automatically, without manual overrides. This host's
`pip --target` layout does not install console scripts, so the run used the
existing `pymacruise` entry-point script with `PYTHONPATH` pointing only at the
fresh target. The resolved package, NumPy, SciPy, nibabel, numba, and ANTsPy
imports came from that target.

The synthetic BASE run completed all eight checked stages, including the real
ANTsPy N4 stage. It wrote 24 final outputs; all manifest hashes verified, and
all 13 NIfTIs were finite and RAS. A same-arguments resume verified and skipped
all eight stages. Both TOADS passes used a two-iteration cap and did not reach
natural convergence. This validates dependency resolution, native imports,
packaged atlases, the installed entry point, workflow execution, and restart;
it is not source-equivalence or natural-convergence evidence.

## Reproduction

From the project root, resolve and install all runtime dependencies plus the
ANTsPy extra into a fresh target:

```sh
.venv/bin/python -m pip install --no-build-isolation \
  --target build/clean-install/with-ants-constrained-site '.[ants]'
```

Then set `PYTHONPATH` to that target and run `pymacruise base-run` using the
fixtures under `tests/data/synthetic/toads_registered_first_iteration_cycle/`
and `tests/data/synthetic/base_workflow/maseg.nii.gz`. The tested run used one
registration level/iteration, two TOADS iterations, one Dura iteration, one
ACE diffusion iteration, and one TGDM iteration; its exact settings and
outputs are saved in the ignored run manifest at
`build/clean-install/base-run-constrained/run.json`. Repeat the same command
with `--resume` to verify all eight stages.

The optional-extra bounds are supported by this current host's package set and
the successful workflow run. Cross-platform dependency resolution, including
other operating systems and Python versions, remains unverified.
