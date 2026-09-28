# TOADS atlas spacing source audit

**Date:** 2026-09-28  
**Roadmap:** P18ac complete  
**Validated implementation commit:** `32e6d804e56b312da04d02ec658abac57f36493a`

## Finding

The selected TOADS Java loader uses the text descriptor's dimensions and
resolutions as the algorithm grid geometry. It reads the referenced topology
and shape data as headerless payloads and does not open or consult an image
header/affine. Therefore the 0.9375 mm root-descriptor resolution is the
resolution used by the source algorithm. The topology image's 1.0 mm NIfTI
zoom reflects the converted source header and must remain visible as a
metadata discrepancy; it must not silently override descriptor spacing or
cause an unvalidated resampling.

## Source evidence

The inspected source is
`toads-cruise/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java`,
SHA-256 `df2f2ace677f2222e94c946730adf485d4401a76f7a9d2b302c9e5c8b652596c`.
At lines 963–990, `loadAtlas` parses topology dimensions and resolutions from
the descriptor, stores the resolutions in `rtx/rty/rtz`, and calls
`loadTemplateImage`. The shape block parses each class's dimensions and
resolutions at lines 1049–1077 before calling `loadShapeImage`. The loaders at
lines 818–859 allocate arrays from the supplied dimensions and decode byte or
little-endian float values from the file; they do not read a header. This is
direct source-body evidence, not a runtime execution of the full Java atlas
class.

The atlas loader in `src/pymacruise/toads_atlas.py` keeps the descriptor
spacing and NIfTI affine zooms as separate metadata. Its preflight warning now
states the source rule. `tests/test_atlas_assets.py` checks that the topology
NIfTI retains 1.0 mm zooms, the descriptor retains 0.9375 mm spacing, every
shape descriptor uses 0.9375 mm spacing, and the discrepancy is reported.
The separately recorded payload audit in `docs/toads-port-plan.md` shows the
RAS NIfTIs restore exactly to the 12 selected raw source arrays in source index
order, with no interpolation.

This resolves which spacing metadata drives the selected Java code. It does
not validate physical atlas-to-input alignment, rigid registration, or
segmentation; those remain open before an atlas execution claim.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **268 passed in 50.27 seconds**
- Pytest: 9.1.1; Python 3.13.2
- Tested commit: `32e6d804e56b312da04d02ec658abac57f36493a`
- Run date: 2026-09-28
