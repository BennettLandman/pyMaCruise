# Synthetic base workflow MaSeg input

`maseg.nii.gz` is an intentionally sparse, deterministic integer label map on
the compact 13³ TOADS loop phantom's T1 grid. It contains all eight selected
MaCRUISE brainstem/cerebellum cleanup labels. Regenerate it with:

```bash
python tests/data/synthetic/base_workflow/generate_maseg.py
```

This tiny input exercises MaSeg grid validation, the three-map cleanup branch,
and the final surface-guided segmentation-refinement output in a runnable
synthetic workflow. Its sparse anatomy is not an anatomical accuracy test.

## Runnable synthetic base run

From the repository root, run:

```bash
./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/toads_registered_first_iteration_cycle/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_registered_first_iteration_cycle/brain-mask.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_registered_field_refit/atlas-11.txt \
  --maseg tests/data/synthetic/base_workflow/maseg.nii.gz \
  --output-dir /tmp/pymacruise-base-synthetic \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 \
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

Rerun the same command with `--resume` to reuse only stage outputs whose
recorded hashes and input signatures still verify. Outputs include a
`run.json`, per-stage logs/manifests, N4 image, both TOADS classifications and
tissue maps, Dura products, ACE/GVF/cropped inputs, native-grid inner/central/
outer level sets, and final `SEG_refine.nii.gz` plus a rule summary. The test
`tests/test_base_workflow_restartable.py` exercises a clean run, no-op resume,
and recovery after damaging a CRUISE output. The compact TOADS fixture reaches
its iteration cap before meeting the energy threshold; a completed run is not
convergence or Java-equivalence evidence.
