# In-vivo current-seed component audit — 2026-10-11

This is a read-only audit of saved files from the corrected T1-grid fused base
checkpoint. It did not rerun MaCRUISE, TGDM, or topology correction.

## Inputs and method

All paths below are relative to the `pymacruise` repository root. Binary masks
were labeled with `scipy.ndimage.label` using 3D 6-, 18-, and 26-neighbor
connectivity. Membership-like images use `value >= 126.9`, matching the
selected TGDM initial isovalue. Signed level sets use `value <= 0` for their
inside mask. “Small” means a connected component with fewer than 100 voxels.

| Saved image | Binary-mask rule | 6-connected components (small) | 18-connected components (small) | 26-connected components (small) |
|---|---|---:|---:|---:|
| `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/wmfill_cropped.nii.gz` | `>= 126.9` | 216 (215) | 62 (61) | 51 (50) |
| `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/wmfill_topology_corrected.nii.gz` | `>= 126.9` | 155 (154) | 4 (3) | 2 (1) |
| `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/inner.nii.gz` | `<= 0` | 155 (154) | 40 (39) | 24 (23) |
| `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/central.nii.gz` | `<= 0` | 155 (154) | 7 (6) | 6 (5) |
| `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/outer.nii.gz` | `<= 0` | 155 (154) | 1 (0) | 1 (0) |

The current crop is 141×173×131 voxels. Thresholded topology-corrected WM fill
has 155 six-connected components, 154 smaller than 100 voxels; the largest is
408,509 voxels. The saved pre-correction WM fill has 216 components, 215
smaller than 100; the largest is 403,338 voxels. Thus the saved Python
topology-correction output reduces the six-connected specks but leaves many.
This audit alone does not establish whether those remaining islands match the
selected Java `AlgorithmTopologyCorrection` output.

## Held-fixed Java TGDM masks

The Python full-grid level sets were cropped to the inclusive provenance bounds
`x=30..170, y=58..230, z=105..235` before comparing them with these saved Java
outputs:

- `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/java-maseg-tgdm-replay/output/inner.nii.gz`
- `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/java-maseg-tgdm-replay/output/central.nii.gz`
- `build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/java-maseg-tgdm-replay/output/outer.nii.gz`

With the same `<= 0` sign rule, Java and Python inner and central masks are
voxel-identical (0 differing voxels). Outer masks differ at 53 voxels. Each
Java mask has the same six-connected component and small-component counts as
its corresponding Python mask: 155 total, 154 under 100 voxels. This supports
the localization of the small islands to the TGDM seed/input, rather than a
Python-versus-Java TGDM discrepancy. It does not show whether the upstream
seed itself matches Java topology correction.

## Surface and label observations

The saved quality audit at
`build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/human-quality-audit-with-refined.json`
reports 158 mesh connected components for each surface. The inner mesh has
157 non-largest components totaling 954 vertices; nearest-voxel sampling maps
942 of those vertices inside the supplied brain mask. The sampled MaSeg labels
are mostly 44 (386 vertices) and 45 (504 vertices). This is a nearest-sampled
vertex diagnostic, not a voxelwise component-to-label assignment. Mesh
connectivity can differ from connected components of a thresholded voxel mask,
so the 155 volume components and 158 mesh components are not directly
interchangeable counts.

The exact label and mask inputs to that audit are
`build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz`
and
`build/invivo-base-2026-10-10-inputs/ras-grid/T1_bet_mask_ras.nii.gz`.

The same audit reports raw-MaSeg terminal support of 81.8786% for GM and
85.8229% for WM, versus 89.2996% GM and 96.7904% WM after source-faithful
segmentation refinement. These are containment/support fractions, not Dice
scores. The findings identify inherited seed fragmentation and a separate
surface-versus-MaSeg support gap; they do not establish the cause of that gap.

## Minimal reproduction

From the `pymacruise` repository root, this reads only the five saved NIfTIs
above and prints thresholded component counts. It does not run a workflow:

```sh
PYTHONPATH=src .venv/bin/python - <<'PY'
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage

root = Path("build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11")
cases = {
    "wmfill_cropped.nii.gz": ("wmfill_cropped.nii.gz", 126.9),
    "wmfill_topology_corrected.nii.gz": ("wmfill_topology_corrected.nii.gz", 126.9),
    "inner.nii.gz": ("inner.nii.gz", 0.0),
    "central.nii.gz": ("central.nii.gz", 0.0),
    "outer.nii.gz": ("outer.nii.gz", 0.0),
}
structures = {
    6: ndimage.generate_binary_structure(3, 1),
    18: ndimage.generate_binary_structure(3, 2),
    26: np.ones((3, 3, 3), dtype=bool),
}
for label, (filename, threshold) in cases.items():
    values = np.asarray(nib.load(str(root / filename)).dataobj)
    mask = values >= threshold if threshold else values <= threshold
    counts = []
    for connectivity, structure in structures.items():
        component_labels, count = ndimage.label(mask, structure=structure)
        sizes = np.bincount(component_labels.ravel())[1:]
        small = int(np.count_nonzero(sizes < 100))
        counts.append(f"{connectivity}-connected={count} (under100={small})")
    print(f"{label}: " + ", ".join(counts))
PY
```

The saved Java/Python mask comparison is recorded in
`build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/java-maseg-tgdm-replay/comparison.json`;
its reproduction command and metric definitions are in
[`in-vivo-fused-input-tgdm-held-fixed-java-comparison-2026-10-11.md`](in-vivo-fused-input-tgdm-held-fixed-java-comparison-2026-10-11.md).
For a direct mask check from saved files, use the following crop and count:

```sh
PYTHONPATH=src .venv/bin/python - <<'PY'
from pathlib import Path
import nibabel as nib
import numpy as np

root = Path("build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11")
crop = tuple(slice(lo, hi + 1) for lo, hi in ((30, 170), (58, 230), (105, 235)))
for name in ("central", "inner", "outer"):
    python_mask = np.asarray(nib.load(str(root / f"{name}.nii.gz")).dataobj)[crop] <= 0
    java_mask = np.asarray(nib.load(str(root / "java-maseg-tgdm-replay" / "output" / f"{name}.nii.gz")).dataobj) <= 0
    print(f"{name}: {np.count_nonzero(python_mask != java_mask)} differing voxels")
PY
```
