# Superseded atlas-aligned synthetic PLUS loop report — 2026-10-12

**Do not use the loop-quality values in this report as current validation.** The fixture's topology-to-MaSeg mapping incorrectly grouped atlas structures 14, 16–18, and 22 as cortical GM and structure 24 as cortical WM. Source-derived MaSeg grouping maps topology class 15 alone to cortical GM 100, classes 25–27 to cortical WM 44, classes 5–6 to CSF 5, topology 1 to background 0, and the remaining listed structures to other tissue 75.

The generator, committed NIfTI label map, provenance, and regression test now use the corrected mapping. All iteration 2–12 ribbon-support results previously recorded for this fixture are superseded. The Numba sampler and selected Java method comparison remain useful method-level evidence, but the large-grid comparison used registration state created from the old labels; it does not validate the corrected full workflow.

See [the corrected source-label terminal run](plus-atlas-consistent-corrected-terminal-2026-10-12.md) for the corrected full nine-stage synthetic run, natural convergence at iteration 10, recursive manifest/hash and RAS checks, explicit ribbon audit, and same-levelset JIST-versus-Lewiner comparison.
