# Synthetic second-pass TOADS preparation for CRUISE (P18cybt)

## Ported behavior

The selected base route after second-pass TOADS can now be called as one
Python API/CLI stage. It optionally removes MaSeg brainstem/cerebellum labels
from GM, WM-fill, and WM-mask; derives and applies one GM-support crop; then
topology-corrects cropped WM-fill using cropped WM-mask as paint. The default
base topology parameters are (6,18), threshold 0.5, levels [0,1], and
regularization 1e-4. Multi-atlas refinement is excluded, matching the selected
base wrapper default. Four RAS NIfTIs are written, including the corrected WM
field consumed by ACE/CRUISE preparation.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_cruise_preparation.py tests/test_toads_tissue_crop.py tests/test_toads_cerebellum.py tests/test_topology.py
```

Result: **19 passed, 0 failed, 0 skipped in 1.29 seconds**.

The committed synthetic fixture verifies optional cleanup of all eight
selected labels, source stage order, shared inclusive GM bounds, cropped RAS
geometry, topology output type, CLI/API array equality, and provenance. The
same constituent cropper and topology kernels have separate selected-Java
reference checks. This composition itself is not a combined Java wrapper run.

Arbitrary JIST orientation-based MaSeg transfer, optional multi-atlas
refinement, ACE/GVF/TGDM execution, and complete base-workflow equivalence
remain open. P18cybt remains Started.
