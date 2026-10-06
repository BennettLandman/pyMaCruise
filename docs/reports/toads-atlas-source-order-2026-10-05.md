# TOADS atlas source-order loading — 2026-10-05

## Implemented behavior

The selected VandyCRUISE atlas images are stored as RAS-oriented `.nii.gz`
files. Their conversion audit records the original axis codes and dimensions.
`load_toads_atlas_data` uses that audit to restore the original Java voxel
index order with axis permutation/flips only; it does not resample values.
For standalone NIfTI atlases without a conversion audit, it preserves the
image's native array order. `toads_alignment_startup` uses the audited loader
for the bundled source atlas. Public subject inputs and outputs remain RAS.

## Validation

Focused command from the PyMaCRUISE repository root:

```bash
./.venv/bin/python -m pytest -q --tb=line tests/test_atlas_assets.py tests/test_toads_alignment_startup.py
```

Result: **14 passed, 0 failed, 0 skipped in 5.61 seconds**. The loader tests
check restoration of asymmetric source-order voxels for audited and
unaudited atlases; the alignment-startup tests check that the bundled atlas
data enters the selected TOADS setup. Existing selected Java comparisons
cover crop bounds, transformed topology, objective initialization, and the
first registration transform; the ordering-carry report separately records
the subsequent growth, closest-label, thinning, and iteration-state
comparisons.

The full suite for the implementation passed on commit
`c8a33033e823d9a2fbce5170ee7ebd984008e801`: **1,202 passed, 0 failed, 0
skipped, 11 expected RAS warnings in 1,187.29 seconds**, peak sampled RSS
**13.539 GB** under the 24 GB watchdog. See
[full-suite report](full-pytest-2026-10-05-toads-ordering-carry.md).

The latest brain-resolution synthetic base workflow, which uses the bundled
atlas through the selected TOADS startup, completed all six stages with both
TOADS passes converged and resumed by verifying/skipping all stages. Its
peak sampled RSS was 6.887 GB under the 16 GB watchdog. See the
[base workflow report](base-workflow-ordering-carry-2026-10-05.md).

## Evidence boundary

This closes the atlas payload loading and source-order mapping task. It does
not establish whole-TOADS or whole-MaCRUISE Java/MATLAB equivalence, nor
real-MRI validation.
