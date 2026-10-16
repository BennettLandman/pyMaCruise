# Exact-input human ACE comparison against selected Java source — 2026-10-16

## Scope and result

The selected MaCRUISE Java ACE harness and Python's public `ace_enhance`
entry point were run on the exact GM, WM-fill, and CSF-boundary arrays used by
the human BASE workflow's CRUISE core. They shared a `(140, 169, 134)` RAS
crop and identical values. The Java wrapper settings were 25 diffusion
iterations, CSF merge placement 2 (after diffusion), GM intensity threshold
126.9, tensor alpha 1.6f, the initial guide Gaussian alpha 1.8f, and the
configured speed Gaussian alpha 1.8.

All delivered ACE products match bit-for-bit:

| Product | Changed voxels | Maximum absolute error | Foreground Dice |
| --- | ---: | ---: | ---: |
| Enhanced GM | 0 | 0 | — |
| Cleaned CSF skeleton | 0 | 0 | 1.0 |
| Thinned CSF skeleton | 0 | 0 | 1.0 |

The unchanged pre-regeneration core-stage products also match these corrected
Java outputs exactly. An initial comparison attempt had re-cropped the
workflow's already-cropped CSF boundary; that comparison is discarded. This
report uses the saved production boundary on its original crop grid.

This is held-fixed ACE equivalence on human-derived inputs. Java and Python
share the upstream arrays created by the Python workflow, so this is not
whole-workflow source equivalence or human anatomical acceptance.

## Input provenance

All inputs and outputs are local to the ignored run directory
`build/sample-singularity-base-source-profile-2026-10-15/java-ace/`.

| Input | SHA-256 |
| --- | --- |
| Cropped GM | `cdb0abec04be4003fe282499a5f4126fcb2a251af5ec521d4c3727ca55382d32` |
| Cropped WM fill | `658c54d5891d6762e3cf19827d2d71e0ee05d61a447a3f7d3d26ba5620ee80ea` |
| Production cropped CSF boundary | `26a5767cd72ffef66587dc1c9f20e877d442fbc53c1e93d3b74d6f1c8e0b121a` |

## Reproduction

Generate selected Java outputs with the source harness:

```sh
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/generate_java_ace_reference.py \
  --place-to-add-csf 2 --diffusion-iterations 25 \
  --input-root build/sample-singularity-base-source-profile-2026-10-15/java-ace/input \
  --output-dir build/sample-singularity-base-source-profile-2026-10-15/java-ace/reference
```

The input directory contains `gray_matter.nii.gz` from the core stage's
`gm_cropped.nii.gz`, `white_matter.nii.gz` from `wmfill_cropped.nii.gz`, and
`multi_atlas_csf.nii.gz` directly from
`base-tissue-preparation/attempt-01/boundary_csf.nii.gz`. The latter is
already cropped by the workflow and must not be cropped again.

Run the public Python API on the same three files:

```python
from pathlib import Path
import nibabel as nib
from pymacruise.ace import ace_enhance

root = Path("build/sample-singularity-base-source-profile-2026-10-15/java-ace")
result = ace_enhance(
    nib.load(root / "input/gray_matter.nii.gz"),
    nib.load(root / "input/white_matter.nii.gz"),
    multi_atlas_csf=nib.load(root / "input/multi_atlas_csf.nii.gz"),
    intensity_threshold=126.9,
    tensor_alpha=1.6,
    gaussian_alpha=1.8,
    diffusion_iterations=25,
    use_updated_csf_in_ace=2,
    merge_method="Merge",
)
```

Compare its three result arrays with `enhanced_gm_java.nii.gz`,
`csf_skeleton_java.nii.gz`, and `thinned_skeleton_java.nii.gz` in the Java
reference directory. The production Python wrapper now explicitly preserves
the selected Java `1.8f` initial-guide and `1.6f` tensor-alpha call
semantics. Focused validation after that change: **22 passed**, peak sampled
RSS **0.333 GB** under an 8 GB watchdog.

The Java output file hashes are enhanced GM
`03fbb7f9193c026b083fdfa382bd5aaf6bf27fea5ae343e6e54f447ff1895d68`, skeleton
`2062ec40aa0039a474ba85b5537cd5fe52299dc00010c769d119e90b5c13af14`, and
thinned skeleton
`c6b7e2a3317ffce0223389c753d7a631911d6621670fe86b7270c997cc73fbd4`.
