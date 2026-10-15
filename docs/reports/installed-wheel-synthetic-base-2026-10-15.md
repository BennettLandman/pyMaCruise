# Installed-wheel synthetic BASE smoke run — 2026-10-15

## Result

Built `pymacruise-0.1.0-py3-none-any.whl` from tested tree commit
`b3b305c5dfbb7473dd3cadfe45a5ce153b7b49f8`, installed it into an isolated
target directory, and ran the installed `pymacruise` entry point with the
working directory outside the source checkout. Python resolved the package
from `build/wheel-smoke/site/pymacruise`, not `src/`.

The synthetic BASE workflow completed all eight stages and recorded 24 final
outputs. A same-arguments `--resume` verified and skipped all eight stages.
Every reported NIfTI is RAS. Peak sampled RSS was 0.732 GB under an 8 GB
watchdog. This is a package/CLI integration smoke test, not a clean dependency
installation or a source-equivalence result: build isolation could not fetch
`setuptools` because network access was unavailable, so the wheel was built
with the already-installed build dependencies, and runtime dependencies came
from the project virtual environment.

Both TOADS passes intentionally used a two-iteration cap. The terminal energy
ratios were `0.0556259789` and `0.118235528`; neither met the selected
convergence threshold. Therefore this run does not claim natural convergence,
anatomical quality, or MaCRUISE equivalence.

## Reproduction

From the repository root, build and install the wheel:

```bash
.venv/bin/python -m build --no-isolation --wheel --outdir build/wheel-smoke/dist
.venv/bin/python -m pip install --no-deps --target build/wheel-smoke/site build/wheel-smoke/dist/pymacruise-0.1.0-py3-none-any.whl
```

The run used the committed fixtures `tests/data/synthetic/toads_registered_first_iteration_cycle/{t1.nii.gz,brain-mask.nii.gz}`, `tests/data/synthetic/toads_registered_field_refit/atlas-11.txt`, and `tests/data/synthetic/base_workflow/maseg.nii.gz`. The exact command is preserved in the ignored run's manifest at `build/wheel-smoke/base-run/run.json`; the installed command was `build/wheel-smoke/site/bin/pymacruise base-run` with one registration level/iteration, two TOADS iterations, one Dura iteration, one ACE diffusion iteration, and one TGDM iteration. Repeating it with `--resume` verified the saved stage outputs.

Run artifacts and the wheel are ignored under `build/wheel-smoke/` and are not committed.
