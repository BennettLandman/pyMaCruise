# Packaged `toads-loop` CLI validation (P18cybn)

## Wheel build and install

Built the standard pure-Python wheel with the already-installed build tools:

```text
./.venv/bin/python -m build --no-isolation --wheel --outdir build/wheel-validation
```

The build used local setuptools 84.0.0. The isolated build first attempted to
fetch its declared `setuptools>=68` dependency, but network access was
unavailable; the no-isolation build completed successfully. The resulting
`pymacruise-0.1.0-py3-none-any.whl` was 134 MB and had SHA-256
`20019572c351589a0ea396a8b26dacf12bc417fba7eb7aef331079ecb92de696`. It was
installed with `pip --no-deps --target` under ignored `build/wheel-validation`
so the wheel itself is not checked into Git.

## CLI validation

The installed entry point reported `pymacruise --help` and
`pymacruise toads-loop --help`. The wheel import resolved to
`build/wheel-validation/site/pymacruise`, and
`get_bundled_atlas_descriptor("cruise-atlas-10obj-toads2008.txt")` found the
packaged default atlas descriptor.

The installed wheel then ran `toads-loop` on the compact committed synthetic
T1/mask/11-class atlas fixture with one registration level, one registration
iteration, and two TOADS iterations. It wrote the native-RAS uint8
classification and four source-scaled Dura inputs. A post-run check verified
all six `.nii.gz` products use the fixture's 13×13×13 grid and RAS geometry.

## Full-suite evidence

The full suite also passed **1,100 tests, 0 failed, 0 skipped, 11 expected RAS
warnings in 2,218.92 seconds** on implementation commit
`a7c142cead48bc13bd1e96854992d903ea3377bf`; see [the full-suite
checkpoint](full-pytest-2026-10-02-n4-cruise-chain.md). The suite run is separate
from the wheel smoke check, which uses the built artifact itself.
