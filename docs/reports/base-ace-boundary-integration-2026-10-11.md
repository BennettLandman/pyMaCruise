# V3 base CSF boundary integration

## Finding

The selected MaCRUISE V3 base wrapper creates a boundary from cropped MaSeg,
removes voxels occupied by second-pass TOADS WM (`>0.01`), then removes voxels
occupied by updated WM (`>25.5`). It passes the result to MaACE with
`UseUpdatedCSFinACE=2`, which merges positive boundary voxels after diffusion
and before Gaussian smoothing. The Python boundary and MaACE operations
already implemented those component behaviors, but the restartable base
routes did not connect them. Thus, the earlier human base run recorded mode 2
without actually applying its boundary map and is not evidence of this part of
the selected V3 workflow.

## Change

`base-tissue-preparation` now saves three cropped RAS maps: raw boundary,
TOADS-cleaned boundary, and the final updated-WM-cleaned boundary. Its
provenance records the source labels, overlooked labels, strict cleanup
thresholds, order, and ACE merge position. The restartable `base-run` and
`base-resume-from-toads` routes both include the final boundary as a hashed
core-stage dependency and pass it to ACE. Core input-grid checks accept only
the matching cropped grid or its full reference grid. Selected base CLI runs that provide preparation provenance reject ACE
mode 2 when the boundary is omitted, and the core API now also rejects any
mode-2 call without a boundary.

## Synthetic validation

Command:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/pytest -q --tb=short \
  tests/test_base_boundary_ace_integration.py \
  tests/test_base_tissue_preparation.py \
  tests/test_toads_cruise_core.py \
  tests/test_base_workflow_restartable.py \
  tests/test_base_resume_workflow.py
```

Result: **19 passed in 26.36 seconds**. Synthetic checks cover boundary
creation and ordered cleanup, exact cropped-grid alignment, handoff into
MaACE mode 2, CLI input provenance, and both restartable workflow routes.
These tests establish synthetic integration only; they do not establish
whole-workflow source equivalence or human-data acceptance. No full-suite
evaluation was run, so no source-line estimate is changed.

## Next check

Rerun the human base tail from its verified second-pass TOADS checkpoint, so
the updated outputs include the missing boundary-to-ACE behavior. Reassess
the ribbon-support metrics and compare against available fixed-component Java
references. The historical Singularity output directory contains no saved
reference products; end-to-end human equivalence therefore remains untested.
