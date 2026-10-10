# Atlas-aligned synthetic plus case

This compact 65×81×62 RAS case downsamples the bundled 2012 LesionTOADS
shape/topology atlas and generates matching T1, FLAIR, brain mask, and MaSeg
labels. Its cortical labels follow the atlas descriptor: topology class 15
(Cerebrum-GM) maps to MaSeg 100 (cerebral GM); classes 25–27 (cerebral WM)
map to 44; CSF classes 5–6 map to 5; and cerebellar, subcortical, and brainstem
classes map to 75 (additional brain tissue). Background is 0. This keeps
noncortical structures out of the cortical GM/WM ribbon checks. Shape maps,
topology, input data, spacing, and descriptor dimensions share the same
physical grid. T1/FLAIR use deterministic tissue intensities; they are
synthetic signals, not clinical data or an independent source-equivalence
reference.

Regenerate the inputs, descriptor, and local resampled atlas maps with
`./.venv/bin/python tests/data/synthetic/plus_vandy_atlas_aligned/generate_fixture.py`.
`provenance.json` records hashes and RAS geometry for each generated NIfTI.
