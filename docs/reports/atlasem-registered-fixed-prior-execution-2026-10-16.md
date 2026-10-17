# AtlasEM registered fixed-prior execution — 2026-10-16

## Implemented handoff

`finish_atlas_em_registered_execution` consumes the descriptor-driven
two-registration result and composes the selected AtlasEM tail:

1. Build the internal mask from strictly positive input channels and MIPAV's
   selected 2.5D hole fill.
2. Transform the registered atlas shapes onto the recropped T1 grid, apply the
   descriptor shape slope, and normalize priors on the internal mask.
3. Compute the selected initial diagonal-deviation estimates, run the initial
   Atlas-prior posterior membership update, then run the fixed-prior
   centroid/covariance/posterior loop with the source membership-change stop
   rule and iteration cap.
4. Generate the hard tissue-index map and the three Dura-input products.

The implementation reuses individually source-checked preprocessing,
registration, prior mapping, EM statistic/update, and output operations. It
does not run the Dura level-set evolution, export NIfTI files, or claim that
these composed calls have been executed against the selected atlas.

## Synthetic validation

The existing 12×12×12 descriptor phantom now runs both rigid registrations
(one iteration each), native-grid recropping, post-recrop membership
initialization, prior transform/normalization, one initial posterior update,
one outer EM iteration, and Dura-input generation. The test verifies class and
grid dimensions, iteration count, and finite Dura products. It is deliberately
capped at one outer iteration; this is synthetic stage-composition evidence,
not natural convergence or optimizer/source workflow equivalence.

The broader focused AtlasEM, registration, Dura, and MIPAV morphology command
passed **43 tests, 3 expected RAS orientation warnings, in 5.98 s** on code
commit `d800db5641f28eefe3677b506c67d06cced61e6a`. The exact command is recorded
in the latest follow-up section of
[`atlasem-focused-regression-2026-10-16.md`](atlasem-focused-regression-2026-10-16.md).

## Remaining evidence gap

The 8-object `ibsr8obj` descriptor and atlas payload are still unavailable.
The Java runtime can execute selected source harnesses but not the full
MIPAV/JIST plugin on this host. Consequently this composition is based on
synthetic descriptor content and is not a selected-atlas or complete
pre-LesionTOADS validation.
