# JIST versus Lewiner folded-cortex surface phantom

This fixture contains three nested, 41³ RAS NIfTI level sets (`inner`,
`central`, `outer`) on a 1 mm isotropic grid. A radial sphere is modulated by
two angular frequencies to create repeated cortex-like folds and ambiguous
cube configurations. The signed level set is sampled on the voxel grid and
stored as float32. Regenerate it with `make_phantom.py`.

The comparison invokes the selected VandyCRUISE/JIST
`IsoSurfaceOnGrid.solveOriginal(image, 18, 0, false)` body and PyMaCRUISE's
scikit-image Lewiner extractor on each identical voxel array. The Java
reference harness uses API stubs only for JIST image/progress and mesh
containers; it compiles and runs the original selected isosurface algorithm
source. It captures the algorithm's index-space point and triangle arrays.
The stubbed normal generator is deliberately outside the measured geometry.

The fixture does not establish whole-workflow equivalence. Vertex nearest
distances depend on tessellation density. Mesh topology and downstream
closest-point thickness summaries are reported separately. The thickness
comparison applies the same PyMaCRUISE closest-point measurement to each
triangulator's three meshes, isolating the impact of extracted geometry; it is
not a direct Java `ClosestPointThicknessForInner` comparison.
