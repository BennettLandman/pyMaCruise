# In-vivo JIST surface extraction comparison — 2026-10-10

## Result

The selected JIST `IsoSurfaceOnGrid.solveOriginal(ImageDataFloat, 18, 0.0f, false)` and PyMaCRUISE source-derived extractor were run on the exact same human inner, central, and outer level sets from the corrected source-profile base run.

| Surface | Faces in both meshes | Connected components in both | Ordered faces identical | Symmetric mean paired-vertex distance | Maximum paired-vertex distance |
|---|---:|---:|---|---:|---:|
| Inner | 444,652 | 209 | Yes | 0.00000212 mm | 0.00000762 mm |
| Central | 466,992 | 209 | Yes | 0.00000209 mm | 0.00000762 mm |
| Outer | 599,544 | 209 | Yes | 0.00000212 mm | 0.00000762 mm |

For all three meshes, vertex mapping is bijective and every face index is identical. The mesh graphs have no boundary or nonmanifold edges. Because corresponding triangles and vertices are matched, the measured maximum paired-vertex displacement bounds the continuous corresponding triangle displacement; it is well below the 1 mm Hausdorff and 0.01 mm mean surface-position thresholds. This is a held-fixed surface-stage comparison only; it does not provide segmentation Dice or compare against historical MaCRUISE surfaces.

The 209 connected components also occur in the selected Java result on the same level sets. This rules out the Python mesh extractor as the cause of those shells. Their cause remains in the input level sets or an upstream stage; the surface comparison does not establish whether every shell is anatomically valid.

## Reproduce

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl .venv/bin/python \
  scripts/compare_workflow_jist_surfaces.py \
  --core-dir build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/cruise-core \
  --work-dir build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/java-reference-surfaces \
  --evidence-class in-vivo
```

The machine-readable results, Java mesh arrays, float32 level sets, and run log remain in the ignored work directory above. Level sets have shape `(208, 309, 305)`, 1 mm spacing, and RAS orientation. Their float32 C-order SHA-256 digests and the JIST source/harness digests are in `comparison.json`. The script's evidence-class flag keeps these results labeled as in-vivo instead of synthetic-only.

This result demonstrates Java-source agreement for surface extraction on one real human input. It is not whole-pipeline equivalence, clinical validation, or evidence that the upstream level sets pass ribbon-quality checks.
