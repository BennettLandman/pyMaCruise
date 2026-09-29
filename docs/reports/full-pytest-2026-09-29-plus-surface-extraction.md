# Full test suite: MaCRUISEplus nested surface extraction

- Date: 2026-09-29
- Tested implementation commit: `e8e95b31b368d0b4bfc22ac0dc8346aceeac0c04`
- Command: `.venv/bin/pytest -q`
- Result: **483 passed in 78.72 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

`extract_plus_surfaces` routes the full-grid inner, central, and outer TGDM
level sets through the existing surface extractor. On the committed synthetic
17³ anisotropic RAS phantom with nested 3, 4, and 5 mm spheres, every mesh
matches direct `extract_levelset_surface` output exactly, has a closed edge
manifold, and reports RAS+ millimetres. Mean physical radius increases from
inner to central to outer.

Focused tests: 8 passed in 0.54 seconds, including the existing surface tests.
The level sets, generator, provenance note, and SHA-256 manifest are committed
in `tests/data/synthetic/plus_surfaces/`. This validates routing and mesh
format/geometry through the existing scikit-image Lewiner substitution. Exact
equivalence to JIST `IsoSurfaceOnGrid` remains a tracked risk and is not claimed.
