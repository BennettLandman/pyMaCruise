# AtlasEM internal mask construction — 2026-10-16

`AtlasEMSegmentation` creates a statistical working mask by requiring every
input channel to be greater than zero, then fills holes with MIPAV's
`AlgorithmMorphology25D.FILL_HOLES`. The Python
`build_atlas_em_internal_mask` follows that order and reuses the already
source-compared MIPAV 2.5D fill operation. This is the AtlasEM internal mask;
the function does not consume or replace the user-supplied MaCRUISE brain mask.

Source hashes:

- `VandyCRUISE/src/edu/jhmi/rad/medic/methods/AtlasEMSegmentation.java`:
  `b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a`
- `mipav-src/src/gov/nih/mipav/model/algorithms/AlgorithmMorphology25D.java`:
  `d3e271151bce7e635777764648ef21c0253485699e577897c4ff4821f75c70f6`

A synthetic two-channel phantom verifies strict all-channel positive support,
hole filling, and the source's cleared in-plane border. The focused preprocessing
and MIPAV morphology command passed **6 tests in 1.38 s** on implementation
commit `0ca74cd69bd2452fc1a4e83f9e1b0c5c77231c7b`. The broader AtlasEM,
registration, Dura, and morphology group passed **43 tests with 3 expected RAS
warnings in 6.50 s** on that commit. The latter group includes the existing
Java-excerpt comparison for the MIPAV hole-fill operation; the new constructor
composition itself is checked on the synthetic phantom.

This adds mask construction only. It does not validate selected atlas
registration, the missing `ibsr8obj` assets, or the full AtlasEM/Dura stage
sequence on human data.
