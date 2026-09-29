# MaCRUISEplus ACE after-Gaussian CSF merge reference

This fixture reuses the committed synthetic 24³ RAS concentric ACE phantom
and runs the selected Java `MaAnatomicallyConsistentEnhance.solve` with
`placeToAddCSF=3` (`AfterGaussian`) and `ReplaceMethod=0` (`Merge`). It stores
the after-Gaussian speed field and downstream ACE outputs. The source first
maps diffused CSF to speed, applies the alpha-1.8 recursive Gaussian, then
sets speed to 0.1 at positive updated-CSF voxels.

The images are synthetic test data only. Regenerate with:

```sh
.venv/bin/python scripts/generate_java_ace_reference.py \
  --place-to-add-csf 3 \
  --output-dir tests/data/synthetic/ace_after_gaussian/java_reference
```

The existing `tests/data/synthetic/ace/java_reference/` fixture remains the
independent mode-2 (`AfterDiffuse`) base-path reference.
