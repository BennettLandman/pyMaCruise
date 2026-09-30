# P18co, P18cp, P39zzz, and P39zzza validation — 2026-09-30

## Full-suite evidence

- Command: `.venv/bin/python -m pytest -q`
- Tested implementation commit: `846232160e3b13b45c4d6cd68b165d1743597ebb`
- Result: **660 passed, 2 expected RAS reorientation warnings, 169.96 seconds**.
- The warnings come from `test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`, which deliberately provides LPS-oriented topology and paint-mask images and verifies conversion to RAS.
- Combined focused validation: `.venv/bin/python -m pytest -q tests/test_toads_dura_outputs.py tests/test_toads_continuing_iteration_slice.py tests/test_lesion_relations.py tests/test_lesion_brainstem_distance.py tests/test_lesion_toads_distance.py` — **25 passed**.

## P18co — TOADS Dura-removal products and wrapper encoding

`generate_toads_dura_removal_outputs` reproduces the selected Java CSF, GM,
white-matter-fill, and white-matter-mask arrays in source class order. The
adapter also restores hard classification into a signed-byte buffer with
background value 1, restores the four float buffers with zero background, and
applies the selected `MedicAlgorithmToads08` wrapper's unconditional float32
multiplication by 255. Raw and wrapper-scaled arrays remain available as
separate stages.

The asymmetric synthetic case matches the Java method excerpts exactly for
classification, all four raw buffers, and all four scaled buffers. It also
preserves the source's CSF alias overwrite order and `CerebellarGM)` spelling
quirk. The selected lower-case-option wrapper is
`toads-cruise/src/.../MedicAlgorithmToads08.java`, SHA-256
`1153b5aa063a3d3f97904494c5c5097c8546c657810d1f3e8044562d24518261`; its
`scaleData` parameter and condition are commented out, so the active ×255
calls are unconditional. The Java harness models the named
`ImageDataMath.scaleFloatValue` operation with direct Java float multiplication
because that dependency implementation is absent from the available source
checkout. It does not exercise the complete MIPAV `ModelImage` lifecycle or
NIfTI serialization. Fixture hashes, selected VandyCRUISE method provenance,
and the wrapper-copy distinction are in
[`tests/data/synthetic/toads_dura_removal_outputs/README.md`](../../tests/data/synthetic/toads_dura_removal_outputs/README.md).

## P18cp — one continuing TOADS iteration slice

`toads_continuing_iteration_slice` composes the source order for one supported
non-stopping case: topology-preserving thinning, an ordinary-OBJECT membership
sweep using the carried pre-iteration closest-label map, the float32
energy-ratio/control decision, populated growing, then relation/boundary/
closest-label refresh. A folded asymmetric synthetic phantom stores four true
RAS-oriented `.nii.gz` inputs, and the test asserts their `('R','A','S')` axis
codes.

The test compares every output field with a deterministic Python integration
snapshot and verifies that the iteration does not stop and that both thinning
and growing process events. This is a wiring/state-handoff regression, not an
integrated Java output reference. The selected component kernels have their
own Java-reference tests; the slice excludes centroid/variance dispatch,
atlas registration, optional inhomogeneity correction, lesion/outlier paths,
and complete initialization/convergence. See the
[`toads_continuing_folded_slice` fixture README](../../tests/data/synthetic/toads_continuing_folded_slice/README.md)
for exact boundaries.

## P39zzz — LesionTOADS relation construction and diffusion

`lesion_build_relation_list` and `lesion_compute_relations` match selected
`LesionToadSegmentation.buildRelationList()` and `computeRelations()` method
excerpts. Six non-cubic synthetic label arrays compare initial/final relation
maps, relation memberships and adjacency rows, relation counts, source pass
counts, unchanged border values, and input immutability. They exercise nested
six-face scan insertion, the source X−/X+/Y−/Y+/Z−/Z+ selection order, and
synchronous diffusion at zero, one, and two limits. Source file SHA-256 is
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`; the
fixture pins its Java excerpt and generated JSON hashes. This does not validate
topology-changing lesion updates or full LesionTOADS execution.

## P39zzza — LesionTOADS brainstem distance input

`lesion_toads_brainstem_distance` builds a seed mask using exact selected
`String.equals` names `Brainstem`, `InterVentricular-WM`, and `FornixWM`, then
reuses the existing Java-matched six-connected capped byte-distance kernel.
The pinned Java method-excerpt comparison covers border and interior seeds,
all exact names, punctuation/case controls, an empty seed mask, and caps 1 and
4 on asymmetric true-RAS `.nii.gz` segmentations. Physical spacing is ignored
by the selected source routine. This validates the name-to-distance adapter
for supplied labels, not LesionTOADS membership evolution or full workflow
equivalence; see the
[`lesion_brainstem_distance` fixture README](../../tests/data/synthetic/lesion_brainstem_distance/README.md).
