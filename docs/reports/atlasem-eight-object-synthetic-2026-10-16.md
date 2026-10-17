# Eight-class AtlasEM synthetic workflow — 2026-10-16

Implementation/test commit: `53d8915b29fd8acc9a214d0e7685e33bac2dca69`

## What was exercised

The restartable `pymacruise atlasem-run` CLI now has a generated eight-class
atlas case in addition to the compact three-class case. The fixture supplies
an eight-class descriptor, topology image, seven positive object priors, T1,
and a matching T1-grid brain mask. It runs both rigid registration stages,
the fixed-prior EM loop for two updates, product restoration, and the
pre-LesionTOADS Dura pass. The test checks the completed stage manifest, the
8-channel native-grid RAS membership product, finite values, and the expected
13 compressed NIfTI products.

All image assets are generated in the pytest temporary directory. This
exercises eight-class orchestration and file contracts; it does not validate
the absent selected `ibsr8obj` descriptor or atlas payload, selected-atlas
registration behavior, or direct Java optimizer parity.

## Validation

Focused AtlasEM CLI and registration regression:

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_cli_workflow.py \
  tests/test_atlas_em_registration_from_descriptor.py --tb=short
```

Result: **6 passed in 1.76 s**.

Full suite on the same implementation commit:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- .venv/bin/pytest -q --tb=short
```

Result: **1,646 passed, 4 skipped, 0 failed, 21 expected orientation warnings
in 728.57 s**. Peak sampled RSS was **14.788 GB** under the 16 GB watchdog; the
run was not terminated for memory use.
