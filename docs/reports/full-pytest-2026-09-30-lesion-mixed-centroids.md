# Full-suite evidence — P39zzy mixed LesionTOADS centroids

- Date: 2026-09-30
- Implementation commit tested: `fbecc45a357576d0c5a583d366d4dc5e1c3d6087`
- Command: `.venv/bin/pytest -q`
- Result: **628 passed**, 2 expected RAS reorientation warnings, in 118.15 seconds.
- Focused validation: 2 mixed-centroid tests passed; the six focused centroid suites
  (mixed, generic, optimized tissue, OUTLIER, ordinary auxiliary, and optimized
  auxiliary) passed 12 tests together.

The committed RAS NIfTI phantom contains all selected class categories in one
13-class vector, including a deliberately absent class, two channels, bias
fields, an OUTLIER map, and both auxiliary maps. For all five update modes,
with and without fields, the Python mixed-category API matches the selected
Java `computeCentroids()` harness and upstream `PowerTable` in raw float32 bits
for class, lesion, and black-hole centroids. Cross-category linked groups,
T1-only optimized-CSF weighting, FLAIR fallback, outlier suppression, ordinary
and optimized auxiliary routing, tie precedence, RAYLEIGH accumulation, and
missing-input validation are covered. The selected Java source reference is
the checkout whose SHA-256 is recorded in the fixture provenance note.

The only warnings are the established tests that deliberately supply LPS
topology volume and mask inputs; both are reoriented to RAS as required.
