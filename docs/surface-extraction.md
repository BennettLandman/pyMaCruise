# Scientific level-set surface extraction

## Production method

MaCRUISEplus calls the VandyCRUISE `GetMeshSurfaceFromLevelSet` wrapper, which calls JIST `IsoSurfaceOnGrid.solveOriginal` with connectivity 18, iso-level zero, and `inclusive=false`. Python production extraction uses a compact 256-case table generated from that supplied Java implementation. It preserves the selected near-zero perturbation, 18/6 case templates, grid-edge interpolation, shared-edge vertex IDs, source traversal order, and face winding. The Java harness is used only for development validation; runtime extraction is Python-native and does not depend on Java or scikit-image.

`pymacruise.surfaces.extract_levelset_surface` normalizes non-RAS input axes with a warning, applies the source perturbation and source-derived templates, and transforms voxel vertices through the complete normalized NIfTI affine. `SurfaceMesh` coordinates and exported legacy VTK `POLYDATA` are RAS+ millimetres. LPS VTK inputs are converted with a warning and reversed face winding; untagged legacy VTK requires an explicit input frame.

## Validation

The selected Java implementation generated all 256 corner-sign patterns used by the Python table. The exhaustive and saved-volume checks show exact ordered triangle connectivity; paired coordinates differ only by floating point rounding (maximum **7.15×10⁻⁶ mm** in the saved base and PLUS surfaces). See [the case-template port report](reports/jist-connectivity18-case-template-port-2026-10-08.md).

The corrected terminal synthetic PLUS case also compares Java JIST and the Python production extractor on identical level-set arrays. Each of the three meshes has identical topology and face connectivity in both results. Both show four connected components and Euler characteristic −8, confirming that those surface properties originate in the level sets for this case. See [the terminal comparison](reports/plus-terminal-jist-surface-comparison-2026-10-12.md).

Synthetic analytic checks cover spheres and ellipsoids, including closed two-manifold edges, Euler characteristic, winding, millimetre bounds, oblique affine handling, RAS normalization, LPS conversion, and VTK round trips. Historical scikit-image comparisons remain useful for documenting the rejected alternative but do not describe the production method.

## Scope

Only ASCII legacy VTK triangular `POLYDATA` is read and written here. This module handles scientific mesh data, not visualization or rendering. The repository uses NumPy/SciPy for mesh geometry and closest-point thickness; it has no VTK/PyVista rendering requirement. Whole-workflow Java equivalence and real-human PLUS acceptance remain separate validation targets.
