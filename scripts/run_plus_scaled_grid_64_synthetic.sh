#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
mkdir -p build/plus-scaled-grid-current build/wheel-plus-natural/mpl-cache

export PYTHONNOUSERSITE=1
export PYTHONPATH=src
export MPLCONFIGDIR="$repo_root/build/wheel-plus-natural/mpl-cache"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

exec .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 \
  --interval-seconds 1 \
  --report-interval-seconds 30 \
  -- .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_scaled_grid/t1.nii.gz \
  --flair tests/data/synthetic/plus_scaled_grid/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_scaled_grid/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_scaled_grid/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt \
  --output-dir build/plus-scaled-grid-current/natural-run-64 \
  --registration-iterations 2 \
  --shrink-factors 1 \
  --smoothing-sigmas 0 \
  --registration-sampling-rate 1 \
  --n4-shrink-factor 1 \
  --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 \
  --crop-padding 1 \
  --init-align-iterations 0 \
  --pyramid-levels 2 \
  --lesion-iterations 1 \
  --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 \
  --finite-zero-variance \
  "$@"
