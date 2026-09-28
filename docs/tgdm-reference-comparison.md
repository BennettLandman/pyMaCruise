# TGDM synthetic Java reference comparison

## Method

`tests/java_reference/` contains a standalone test harness built around the
supplied MaCRUISE `MaGenericTGDM` source, VandyCRUISE profile and topology
sources, and its `DistanceFieldyh` reference implementation. The harness
recreates the adapter call sequence recorded in
[`cruise-tgdm-port.md`](cruise-tgdm-port.md) using the committed 25³
concentric-tissue phantom. The Java outputs are checked-in `.nii.gz` fixtures;
no subject images are used. Run `scripts/generate_java_tgdm_reference.sh` to
recreate them.

The Python comparison runs the same default selected wrapper settings. It
checks the membership signed seed, smoothing, inner seed, tissue-sign field,
central evolution, central shrink, final inner evolution, and outer evolution.
Level-set masks use `value <= 0`. Surface points come from the zero-level
triangulated meshes; the reported symmetric mean and Hausdorff use exact
closest points on opposing triangle surfaces, sampled at each mesh's vertices.
The fixture has 1 mm isotropic spacing and identity RAS affine. Run
`scripts/compare_tgdm_java_reference.py` to regenerate the JSON evidence under
`docs/reports/`.

## Results

| Stage | Voxel Dice | Symmetric mean surface distance (mm) | Hausdorff (mm) | Mean absolute level-set error |
|---|---:|---:|---:|---:|
| Initial signed seed | 1.000000 | 0 | 0 | 6.38e-8 |
| Smoothed initial | 1.000000 | 6.97e-5 | 2.34e-4 | 8.64e-4 |
| Initial inner | 1.000000 | 4.45e-7 | 1.88e-6 | 3.20e-7 |
| Central | 1.000000 | 2.86e-7 | 3.81e-6 | 3.45e-7 |
| Central shrink | 1.000000 | 4.43e-7 | 2.86e-6 | 2.50e-7 |
| Final inner | 1.000000 | 5.51e-8 | 6.42e-7 | 1.44e-7 |
| Outer | 1.000000 | 5.72e-4 | 0.0182 | 3.91e-4 |

Every segmentation Dice meets the project minimum of 0.97; all measured mean
surface distances are below 0.01 mm, and all Hausdorff estimates are below
1 mm. The discrete tissue-sign field matches exactly. All outputs are finite,
use matching RAS geometry, and the nested masks have complete voxel
containment on this phantom.

The JSON report records maximum level-set errors as well as MAE/RMSE, per-stage
Dice, and both directed surface means. The outer level set has the largest
maximum scalar error (0.02991); its segmentation Dice remains 1.0 and its
surface metrics pass the agreed spatial thresholds.

The original `solveSmoothSurface` wrapper fast-marches the already signed
membership seed again before its PDE iterations. This second reinitialization
was missing in the first Python draft; adding it reduced the smoothing surface
mean from 0.0586 mm to 0.00012 mm. A regression test now compares all stages
against the saved Java reference volumes.

## Limits

This means the Python implementation matches the Java harness reference for
the selected voxel TGDM wrapper on this compact synthetic shape. It does not
establish broad workflow equivalence, biological segmentation quality, or
behavior on arbitrary grids. The separate committed 41³ folded-cortex fixture
exercises Python narrow-band evolution and topology invariants but has not been
compared with Java references. Mesh extraction and thickness remain separately
validated stages. P29 performance characterization is recorded in
[`tgdm-synthetic-benchmark.md`](tgdm-synthetic-benchmark.md).
