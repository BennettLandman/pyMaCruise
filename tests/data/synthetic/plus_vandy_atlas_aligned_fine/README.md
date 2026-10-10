# Fine-grid atlas-aligned synthetic plus case

This 81×101×77 case resamples the bundled LesionTOADS 2012 atlas onto an
approximately 2 mm RAS grid. It preserves the descriptor's source shape and
topology maps and derives MaSeg labels from named topology classes: 5/6 → CSF
(5), Cerebrum-GM (15) → cortical GM (100), cerebral WM (25–27) → cortical WM
(44), and cerebellar, subcortical, and brainstem regions → additional tissue
(75). T1 and FLAIR are deterministic tissue-coded signals; they are synthetic
inputs and do not model clinical image contrast or establish source
equivalence.

The generator is:

```sh
./.venv/bin/python tests/data/synthetic/plus_vandy_atlas_aligned/generate_fixture.py \
  --output-dir tests/data/synthetic/plus_vandy_atlas_aligned_fine \
  --shape 81 101 77
```

`provenance.json` records the image hashes, data types, shape, RAS axes, and
affines. Use `pymacruise plus-full-run` with this case to exercise the raw-input
workflow; its run manifest records the selected surface-evolution settings.
The explicit ribbon audit additionally checks WM/GM/CSF support, nesting,
closed connected genus-zero meshes, and production triangle and level-set
thickness against a 5 mm analytic target. Workflow completion and Java
comparisons must be reported separately from passing those anatomical-quality
checks.
