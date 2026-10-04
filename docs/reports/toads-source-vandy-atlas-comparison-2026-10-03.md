# TOADS source-atlas comparison — 2026-10-03

## Scope

This is a synthetic Java-source/Python comparison on the committed
`(36,45,34)` Vandy-atlas-derived T1 fixture. The selected Java algorithm reads
headerless atlas payloads in their original descriptor index order. Python's
first comparison used canonical RAS arrays and diverged at crop/startup. The
source-order diagnostic then reproduced the Java crop and startup state. That
behavior is now used by the production TOADS registration path. This validates
the captured source slice only; it does not establish complete TOADS or
MaCRUISE equivalence, convergence on brain-resolution data, or real-data
performance.

## Atlas conversion and ordering

The conversion audit covers 35 bundled atlas images. The source-payload
round-trip check restored all 32 raw payloads exactly using their source
dimensions, datatype, and MIPAV voxel order. The three Analyze images were
voxel-exact after canonicalizing through their source header affine. Their
source-axis metadata still disagrees with the audit (`L,A,S` in the headers vs
`L,P,S` in the audit), so physical origin/direction claims for those examples
remain unverified. Those SPECTRE images are not used in this Vandy comparison.

The implementation now restores audited NIfTI payloads to Java index order
for TOADS registration. Standalone NIfTI atlases without a conversion audit
retain their native array order, using their affine axis codes. The T1 and
final outputs remain on the RAS input grid; the internal atlas index order
preserves the source algorithm's descriptor-based mapping.

## Matched startup and first-cycle results

Java used the original VandyCRUISE 11-structure raw descriptor and the
instrumented `VandyCoreSyntheticRunExport` harness. The synthetic T1 values
were exported in Fortran order to preserve MIPAV's x-fastest voxel convention.
Python used the committed converted NIfTI atlas through the production loader,
one registration level, one registration update per branch, padding 10, and
relation limit 3. Python's trace requires at least two loop iterations; the
captured first iteration is compared with Java's one-iteration run.

- Java and Python crop bounds match exactly: `[3,35] [3,42] [0,30]`; both use
  the `(53,60,51)` padded grid.
- The initial membership objective matches: Java `171.02353`, Python
  `171.02352905273438`.
- The source final rigid transform and Python post-refresh transform agree
  within `1.7e-7` across the six values (the Java diagnostic prints six
  significant digits).
- Transformed shape priors: mean absolute difference `1.0611e-8`, maximum
  `2.8611e-6`; 4,568 of 1,783,980 values exceed `1e-6`.
- Initial memberships: mean absolute difference `1.0744e-9`, maximum
  `1.9670e-6`; 66 of 1,783,980 values exceed `1e-6`.
- First exact membership update: mean absolute difference `5.6475e-10`,
  maximum `3.7849e-6`; 67 of 1,783,980 values exceed `1e-6`.
- Transformed constructor topology matches all 162,180 voxels exactly.
- First thinning labels match all 162,180 voxels exactly.
- Iteration 1 objective and energy ratio match: Java `50.30009 / 1.0909222`,
  Python `50.3000908 / 1.09092224`.

The reported values are selected Java-reference comparisons on one synthetic
input. They are not a Dice/Hausdorff/mean-surface-distance acceptance result
for the final MaCRUISE workflow. Python's final trace output is uint8 on the
input's `(36,45,34)` RAS grid. No real MRI was used.

## Evidence and remaining work

- Java source log: ignored diagnostic file
  `build/java-reference/source-vandy-atlas-first-loop.log`.
- Python log after changing the production atlas loader: ignored diagnostic
  file `build/java-reference/python-source-vandy-postfix-first-loop.log`.
- Focused regression command after the implementation change:
  `PYTHONPATH=src .venv/bin/pytest -q tests/test_atlas_assets.py tests/test_toads_atlas_startup.py tests/test_toads_registered_first_iteration_cycle.py`
  — **14 passed in 9.64s**.
- The deterministic inputs are committed under
  `tests/data/synthetic/toads_downsampled_vandy_atlas/`; the atlas NIfTIs are
  bundled under `data/atlases/`.
- No full pytest suite was run for this change.

Next, continue the Java/Python comparison through later TOADS updates and
convergence, then rerun the restartable base workflow using source-ordered
atlas inputs. Keep those results separate from this first-cycle match.
