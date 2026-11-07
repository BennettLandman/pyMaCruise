# PLUS atlas selection — owner correction, 2026-10-18

> Historical dated record. Current testing readiness, native/Python atlas preparation and run instructions supersede stale open-work statements below; see [the October 27 readiness report](plus-testing-readiness-2026-10-27.md).

The project owner directs that PLUS use `ibsr11objV3` in place of the previously required `ibsr8obj` AtlasEM resource. This supersedes the atlas-selection requirement in older reports, plans, and author-materials requests. The historical MATLAB `run_EM_jist.m` still names `ibsr8obj`; that is evidence of the checked-out historical configuration, not the intended target after this correction.

## Available assets and implementation status

The supplied source tree contains `VandyCRUISE/src/Atlas/ibsr11objV3/LesionToads-11-obj-atlasV3.txt`. It declares 12 structures including background, a 172×216×164 grid, and 0.9375 mm descriptor spacing. All 13 referenced raw payloads (one topology and 12 shape maps) exist relative to the source `Atlas/` directory. Their sizes match signed-byte topology and little-endian float32 shapes at the declared dimensions. The JSON inventory beside this report records their source hashes; the unreferenced small topology is not part of this descriptor's dependency set.

The atlas is available, so obtaining the old `ibsr8obj` package is no longer the atlas prerequisite. The selected Python AtlasEM adaptation is now implemented in `scripts/adapt_ibsr11objv3_atlas.py`. It writes a Python AtlasEM descriptor plus derived RAS `.nii.gz` payloads under a new or empty output directory; source files are left untouched. The implementation:

- Resolves the descriptor's duplicated `ibsr11objV3/` prefix and tolerates the trailing tab in the Thalamus row.
- Converts raw topology and little-endian float32 shapes using the selected Java loader's x-fastest indexing, then applies the declared MIPAV axis orientations to RAS.
- Writes 13 common-grid NIfTI payloads at 0.9375 mm and an AtlasEM-format descriptor preserving the 12 class names, labels, order, intensity profiles, and lesion models.
- Records source/output hashes, input geometry, orientation transform, and limitations in `adaptation-manifest.json`. The misleading source XML geometry is not used.

Recreate the build-only adapter output from the `pymacruise` repository root with:

```sh
PYTHONPATH=.:src .venv/bin/python scripts/adapt_ibsr11objv3_atlas.py \
  --descriptor ../toads-cruise-scm-2026-09-26/VandyCRUISE/src/Atlas/ibsr11objV3/LesionToads-11-obj-atlasV3.txt \
  --output-dir build/ibsr11objV3-python-ras
```

The adapted descriptor is `build/ibsr11objV3-python-ras/LesionToads-11-obj-atlasV3-atlasem-python.txt`. The actual AtlasEM descriptor parser and registration-input loader accept it: class-ordered shape priors load as `(172, 164, 216, 12)`, topology labels match all 12 declared structures, and atlas spacing remains 0.9375 mm. Focused descriptor, registration-adapter, and conversion tests pass (19 tests); `git diff --check` is clean.

This closes the selected-atlas data/descriptor adaptation prerequisite for the Python AtlasEM stage. It does not run the full AtlasEM/PLUS pipeline on human data or establish numerical equivalence against the Java source. The runtime defaults and source-reference driver remain unchanged. The separate LesionTOADS stage retains its existing descriptor selection; this correction applies only to the AtlasEM resource. Whole PLUS source validation and human convergence remain open.
