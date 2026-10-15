# Clean dependency-install attempt — 2026-10-15

## Result

A separate target installation was attempted from the local PyMaCRUISE source
with build isolation disabled. Metadata generation succeeded, but pip could not
resolve the configured package index hostname while requesting NumPy; all
retries ended with `nodename nor servname provided, or not known`, followed by
`No matching distribution found for numpy<3,>=1.23`. No dependency package was
installed into the target directory.

This is an environment/network blocker, not evidence that the package metadata
or dependency declarations are invalid. It leaves clean dependency resolution
and installation unverified on this host. The installed-wheel synthetic BASE
check used the existing project environment's runtime libraries and should be
interpreted only as a wheel/entry-point/workflow smoke test.

## Reproduction

From the project root, with the current virtual environment:

```sh
.venv/bin/python -m pip install --no-build-isolation --target build/clean-install/site .
```

The failure is due to the host not resolving the configured package index; a
network-enabled runner or a complete local wheelhouse is needed to finish this
check.
