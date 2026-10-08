# Brain-resolution thickness impact of surface extraction

This compares the selected Java thickness-wrapper aggregate over Java IsoSurfaceOnGrid meshes with the production Python thickness API over Lewiner meshes extracted from the same committed level sets. Mesh topology differs, so this report compares aggregate values and does not match vertices one-to-one.

| Metric | Selected Java / JIST meshes | Production Python / Lewiner meshes | Difference |
|---|---:|---:|---:|
| Mean thickness (mm) | 1.649155159 | 1.649050572 | -0.000104587 |
| Sample standard deviation (mm) | 0.498031830 | 0.498030609 | -0.000001221 |
| Nonzero inner vertices | 116,165 | 116,178 | — |

The production Python thickness call took `463.4 s` and peaked at `0.631 GB` RSS. No thickness-specific acceptance threshold has been set; the report quantifies the current surface-extraction substitution rather than declaring it equivalent. The P08 surface-position checks and topology differences are documented in the [current-base surface report](current-base-jist-surface-comparison-2026-10-07.md). Reproduce after generating the selected Java comparison with `./.venv/bin/python scripts/compare_jist_closest_point_thickness_brain_resolution.py` by running `./.venv/bin/python scripts/compare_current_base_thickness_substitution.py`.
