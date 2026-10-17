# Synthetic AtlasEM → pre-LesionTOADS Dura run — 2026-10-16

## Runnable evidence

The new test `tests/test_atlas_em_synthetic_full_handoff.py` runs a
restartable-style deterministic synthetic chain from a generated descriptor
and `.nii.gz` atlas payloads:

1. Load class-ordered atlas priors and descriptor intensity profiles.
2. Initialize memberships, perform the initial scale-2 rigid update, uncrop,
   recompute boundaries, and recrop.
3. Reinitialize memberships, derive class bounds, and run the post-recrop
   single-scale rigid update.
4. Build the AtlasEM internal mask and transformed/normalized priors; perform
   the initial Atlas posterior update and one fixed-prior EM iteration.
5. Generate and restore the membership, segmentation, internal-mask and
   three-channel Dura-input products as compressed NIfTI on the native RAS T1
   grid.
6. Feed the saved-product form of the Dura inputs to the selected Dura code and
   complete one outer pass.

Run the focused end-to-end test:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_synthetic_full_handoff.py
```

Result: **1 passed in 0.78 s**. The broader AtlasEM/registration/Dura/morphology
group passed **44 tests, 3 expected RAS orientation warnings, in 6.35 s** on
test commit `c5abb414fe504a0f280ec2806ad8f32c204b66cb`.

## Evidence limits

The synthetic chain completes, writes finite RAS products, and executes one
Dura pass. The registration and EM stages are deliberately limited to one
iteration each; the test does not establish natural convergence, selected
Java optimizer parity, or anatomical quality. The descriptor and atlas payload
are generated phantoms, not the unavailable selected `ibsr8obj` atlas. This is
synthetic workflow evidence, not selected-atlas, real-human AtlasEM, or full
PLUS source equivalence.
