# AtlasEM descriptor metadata and resource inventory

## Scope

Commits `34028ff9a39d78e7f187c8b03eb446fd3dbabcd7` and `ceeadc8bb0c093e3ddc325656801b011e9f2cb36` add a Python reader for the selected VandyCRUISE `Structure Atlas File` text format. It captures descriptor class order, labels and topology names; topology dimensions and resolutions; shape sample count, dimensions, resolutions and orientation; shape-file references; registered-shape vector; intensity and variance profiles; lesion and black-hole model values; modality weights; optimized factors; and shape slope. It follows Java defaults (class zero unregistered, zero intensity profiles, unit variance/weights/optimized factors, zero lesion models, and shape slope 2) and Java registered-shape override rule (`value == 1`). Relative payload references are resolved under the descriptor directory and path escapes are rejected.

The reader inventories topology and shape payload paths and reports missing files. `require_atlas_em_payloads` provides an explicit pre-execution check. Neither function reads voxel payloads or runs AtlasEM. In particular, this work does not establish that the absent selected `ibsr8obj` assets are available or that the AtlasEM workflow is complete.

## Source basis

Selected source: `DemonToadDeformableAtlas.java`, SHA-256 `96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a`. The implementation follows `loadAtlas` section defaults and parsing for `Structures`, `Topology Atlas`, `Shape Atlas`, `Registered Shapes`, `Intensity Atlas`, `Intensity Variance`, `Lesions`, `Black Holes`, `Slope`, `Modality weights`, and `Optimized factors`. Java's raw topology-byte and little-endian shape-float readers remain unported in this slice. The descriptor's per-shape repeated dimensions/resolutions are not used by the Java loader to allocate shape arrays; it uses the topology atlas dimensions.

## Validation

Command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/pytest -q tests/test_atlas_em_descriptor.py
```

Result on implementation commit `ceeadc8bb0c093e3ddc325656801b011e9f2cb36`: **6 passed in 0.68 s**. Synthetic temporary descriptors verify class ordering, spatial metadata, shape and topology path inventory, missing-payload reporting, successful preflight after resources appear, registered-shape defaults and Java override semantics, intensity/profile aliases and defaults, malformed records, and traversal rejection. This is parser behavior validation against source inspection, not Java output comparison, atlas payload validation, or an AtlasEM execution test. No full suite was run.

## Remaining AtlasEM work

The `ibsr8obj` descriptor and payload remain unavailable. Once obtained, bundle them in the repository, confirm all references with the resource preflight, convert any image payload used as a stored project image to `.nii.gz` while preserving provenance, and compare Java-compatible payload loading and selected initialization/registration behavior against the original.
