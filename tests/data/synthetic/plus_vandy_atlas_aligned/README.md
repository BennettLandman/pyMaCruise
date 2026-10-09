# Atlas-aligned synthetic plus case

This compact 65×81×62 RAS case downsamples the bundled 2012 LesionTOADS
shape/topology atlas and generates matching T1, FLAIR, mask, and explicit
MaSeg-like labels (0 background, 5 CSF, 44 WM, 100 GM). Shape maps, topology,
input data, spacing, and descriptor dimensions share the same physical grid.
T1/FLAIR use deterministic tissue intensities; they are synthetic signals,
not clinical data or an independent source-equivalence reference.

Regenerate the inputs, descriptor, and local resampled atlas maps with
`./.venv/bin/python tests/data/synthetic/plus_vandy_atlas_aligned/generate_fixture.py`.
`provenance.json` records hashes and RAS geometry for each generated NIfTI.
