# Lewiner Marching Cubes paper-to-code notes

## Published method and implementation reference

Lewiner T, Lopes H, Vieira AW, Tavares G. “Efficient Implementation of
Marching Cubes’ Cases with Topological Guarantees.” *Journal of Graphics
Tools*. 2003;8(2):1–15. DOI
[10.1080/10867651.2003.10487582](https://doi.org/10.1080/10867651.2003.10487582);
[author-hosted paper and source](https://thomas.lewiner.org/pub/marching_cubes_jgt.html).

The paper gives an implementation of the Marching Cubes 33 cases that resolves
ambiguous configurations and provides topological guarantees for the sampled
isosurface. The scikit-image API identifies its default `method="lewiner"` as
a Cython port of Lewiner's C++ implementation and states its ambiguity
resolution/topological guarantee ([`skimage.measure.marching_cubes`
documentation](https://scikit-image.org/docs/stable/api/skimage.measure.html)).
For the library as a whole, see van der Walt et al., “scikit-image: image
processing in Python,” *PeerJ* 2:e453 (2014), DOI
[10.7717/peerj.453](https://doi.org/10.7717/peerj.453).

## Code crosswalk

| Python code | Relationship to the paper/library | MaCRUISE/JIST boundary |
| --- | --- | --- |
| [`surfaces.py`](../../src/pymacruise/surfaces.py): `extract_levelset_surface`, `marching_cubes(...)`; [`plus_pipeline.py`](../../src/pymacruise/plus_pipeline.py): `extract_plus_surfaces`; [`cli.py`](../../src/pymacruise/cli.py): `surface-extract` command | Calls the external scikit-image Lewiner triangulator (`method="lewiner"`, `step_size=1`, `allow_degenerate=false`) on the shifted zero-level set; plus and CLI wrappers route level sets through this operation. The MC33 case tables and ambiguity handling are not reimplemented in PyMaCRUISE. | The selected MaCRUISE Java wrapper sets connectivity 18 and calls `IsoSurfaceOnGrid.solveOriginal(..., 18, 0, false)`. scikit-image exposes no matching digital-connectivity argument. Lewiner MC33 is therefore a different triangulation path; no JIST-equivalent topology or vertices are claimed. |
| [`surfaces.py`](../../src/pymacruise/surfaces.py): `shifted` near-zero handling | Preserves the selected Java `GetMeshSurfaceFromLevelSet.solveOriginal` behavior: after subtracting the requested level, values in `(-0.1, 0]` become `-0.1` and values in `(0, 0.1)` become `0.1`. | This threshold perturbation is inherited from the Java source contract; it is not part of Lewiner et al.'s published algorithm. Exact zero consequently falls on the negative side. |
| [`surfaces.py`](../../src/pymacruise/surfaces.py): `index_vertices @ affine...` and winding correction | Transforms scikit-image's array-index coordinates with the full normalized NIfTI affine so `SurfaceMesh` vertices are RAS+ millimetres, including oblique orientation and spacing; corrects winding to retain outward orientation. | Physical coordinate conversion and winding normalization are PyMaCRUISE output conventions, not claims made by the Marching Cubes paper. |

## JIST substitution evidence and remaining risk

The selected Java source contract is in `toads-cruise/src/edu/vanderbilt/masi/plugins/CRUISE/utilities/GetMeshSurfaceFromLevelSet.java` and `toads-cruise/src/edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java`: wrapper default connectivity is 18; `solveOriginal` receives iso-level 0 and `inclusive=false`. The Java near-zero perturbation is independently retained by `surfaces.py`. Direct selected-JIST comparisons run both algorithms on identical synthetic level sets:

- Three nested folded-cortex level sets had matching vertex/face/Euler and edge-count invariants. After vertex mapping, only 69.2–70.6% of unordered faces matched. Sampled symmetric vertex-to-triangle means were `1.73e-7`–`1.83e-7 mm`, with maxima `0.83e-6`–`1.86e-6 mm`.
- Adversarial torus, narrow-neck dumbbell, near-touching components, and boundary-crossing sphere cases matched expected component/Euler counts with no nonmanifold edges. Mapped-face agreement was 67.8–71.8%; sampled means were `1.62e-7`–`1.88e-7 mm`, maxima `0.76e-6`–`1.63e-6 mm`.
- On the nested folded surfaces, applying the same Python closest-point thickness operator to both meshes gave a mean thickness difference of `-0.000546 mm`. This isolates geometry's contribution under that Python operator; it is not a comparison against Java thickness code.

The tiny sampled distances do not bound continuous triangle-to-triangle Hausdorff distance, and differing faces establish that local triangulation is not identical even where counts agree. P30 analytic fixtures validate the selected Lewiner path's geometric behavior, not equivalence to JIST. P37 must retain mesh-connectivity and downstream thickness as an explicit substitution risk. See the folded and adversarial source-comparison reports [`surface-extraction-java-comparison-2026-09-30.md`](../reports/surface-extraction-java-comparison-2026-09-30.md) and [`surface-extraction-adversarial-java-comparison-2026-09-30.md`](../reports/surface-extraction-adversarial-java-comparison-2026-09-30.md), plus [`surface-extraction.md`](../surface-extraction.md).
