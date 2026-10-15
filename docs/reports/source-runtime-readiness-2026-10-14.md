# Original full-workflow runtime readiness

Checked from the local PyMaCRUISE host on 2026-10-14 while preparing the fresh-reference work in P36.

## Available

- OpenJDK `21.0.11` and `javac 21.0.11` are on `PATH`.
- The supplied MIPAV and JIST source trees and their dependency JARs are present.
- There are no compiled `.class` files in the supplied MIPAV source tree. The selected Java comparison harnesses compile from source and run locally; this is how the held-fixed TGDM result was produced.

## Not available on this host

`singularity`, `apptainer`, `matlab`, `docker`, `podman`, `ant`, `mvn`, and `gradle` are not on `PATH`; no matching MATLAB or MIPAV application bundle or Singularity image was found in the inspected application and project locations. Therefore the original Singularity-wrapped MATLAB/MIPAV/JIST BASE or PLUS workflow cannot be launched locally from the current installation. This does not block Python workflow execution or selected-source Java harness validation.

The downloaded human output remains a secondary historical reference. The available sidecar v3.2.0 log is for a different `(176, 240, 256)` `sub-010260` input, while the downloaded root output is `(184, 256, 256)`. No fresh complete original workflow output with a matched input/configuration manifest is presently available. P36 remains open for controlled source execution and capture.

## Reproduction

From this host, check executable availability with:

```bash
for exe in singularity apptainer matlab docker podman ant mvn gradle java javac; do
  if command -v "$exe" >/dev/null 2>&1; then printf '%s: %s\n' "$exe" "$(command -v "$exe")";
  else printf '%s: unavailable\n' "$exe"; fi
done
```

The source checkout/JAR inventory used read-only filesystem searches. No files outside `pymacruise` were modified.
