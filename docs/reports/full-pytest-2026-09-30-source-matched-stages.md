# Full test suite: source-matched stages, packaging, and paper annotations

## Run record

- **Implementation commit:** `aa908bac7fbb06eeb31b8a4abb21f7ffad810c3e`
- **Command:** `.venv/bin/pytest -q`
- **Result:** **918 passed, 4 warnings, 618.21 seconds (10:18)**
- **Failures/skips:** none reported

The four warnings were expected orientation notices under the repository's RAS contract. Two were emitted while the calculator CLI test constructed its expected result from LAS-oriented synthetic inputs; two were emitted by the Dura topology test for LPS-oriented topology and paint-mask inputs. The CLI and topology adapters reorient those inputs to RAS and report the conversion.

This is one complete pytest run against the implementation commit above. It includes the source-matched TOADS exact-pass and loop-control slices, LesionTOADS initial-membership-to-thinning composition and numeric fixes, calculator CLI contract, isolated atlas wheel install, surface extraction comparisons, and the code-comment/documentation updates. It does not establish complete TOADS/LesionTOADS workflow equivalence, real-data readiness, or P37 acceptance.

## Focused synthetic evidence within the run

- TOADS selected exact initialization → first thinning → centroids → variances → membership update: **4 passed**; all tested float32 state and objective bits matched the selected Java harness on the folded anisotropic RAS fixture. The fixed-input boundary is described in P18cyf–P18cyi.
- TOADS stop/continue control transition: **7 passed**; the test checks the selected branch ordering and final-growth behavior with registration held fixed (P18cyj).
- LesionTOADS selected initial standard-13 membership and first thinning: **2 passed**; the 17³ folded RAS composition matched all pinned maps, the 1,917-event heap trace, topology decisions, and final state (P39zzzzzzak).
- Affected LesionTOADS numeric/source references after the `PowerTable` and float32 operation-order corrections: **10 passed**. The selected Java table is used with effective float32 exponent 1.2; this evidence is method-slice scoped.
- Generic calculator CLI: **2 passed**, including RAS reorientation warning, `.nii.gz` enforcement, API result parity, and provenance fields (P42c).
- Atlas packaging: **1 passed**, building an sdist and wheel, installing the sdist-built wheel into an isolated target, comparing installed descriptors and every referenced asset SHA-256 with repository sources, checking representative NIfTI orientation, and resolving default, TOADS2012, and LesionTOADS assets without checkout paths (P42b).
- Surface extraction/thickness: **16 passed** across the direct JIST comparison-adjacent tests; the separate comparison reports retain sampled-distance and connectivity limits (P57).
- GVF paper mapping: **8 passed**; LesionTOADS paper-related modules: **22 passed** before the final arithmetic correction and the affected final source-reference subset above passed after that correction. The complete run covers all of them.

The associated paper-to-code evidence is in [paper-annotations-2026-09-30.md](paper-annotations-2026-09-30.md). The surface comparison reports are [folded fixtures](surface-extraction-java-comparison-2026-09-30.md) and [adversarial fixtures](surface-extraction-adversarial-java-comparison-2026-09-30.md).
