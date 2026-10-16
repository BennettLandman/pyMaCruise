# Exact-input human Dura comparison against selected Java source — 2026-10-16

## Scope

This is a held-fixed source comparison of the MaCRUISE V3 Dura stage on the
exact first-pass arrays from the restartable human BASE case documented in
[`in-vivo-root-sample-base-run-2026-10-16.md`](in-vivo-root-sample-base-run-2026-10-16.md).
Both implementations received identical 184×256×256, 1 mm RAS inputs. The
selected Java wrapper was run with five TGDM outer iterations, membership
level 0.3, and GM topology regularization 0.02 (WM 0.05). It is not an
independent full-workflow source execution.

The input image hashes are:

| Input | SHA-256 |
| --- | --- |
| N4 T1 | `c5f56de03e8bbeae23fa1b860abb3772949d44387426d05670202f5041ed160c` |
| WM fill | `4d971c796c508da01e95777506608b84d2ba5ab219c2ca75df62a76e763f0382` |
| GM | `c9eded75cb92631aac25e6ae3993e4f693aaa4baef91cbb7054fd70e469df0f9` |
| WM mask | `1b427a66800a9d88892e11b3261e12f2d19f4f8b95c232dbe680ffa9fed4ab33` |

## Results

The topology-corrected WM, WM+GM union, GM, stripped T1, corrected WM, and
corrected GM are bit-identical across Java and Python. The final brain mask is
also bit-identical (Dice 1.0).

| Product | Maximum absolute error | Changed values | Sign changes / Dice |
| --- | ---: | ---: | --- |
| Inner distance | `9.54e-7` | 166 | 0 / nonpositive Dice 1.0 |
| Outer distance | `9.54e-7` | 161 | 0 / nonpositive Dice 1.0 |
| Dura level set | `2.951e-6` | 94,433 | 0 / nonpositive Dice 1.0 |
| TGDM passes 1–5 | `0.760e-6` to `2.713e-6` | 84,153–94,433 per pass | 0 / nonpositive Dice 1.0 |
| Final brain mask | 0 | 0 | foreground Dice 1.0 |

Floating-point differences in the distance and level-set fields are small and
never change their signs. Thus the Dura stage's output segmentation and mask
match exactly on this human-derived held-fixed input; its continuous fields
are numerically very close but not bit-identical.

## Reproduction

From the repository root, create the Java input pack and run the selected
source harness:

```sh
./.venv/bin/python scripts/pack_dura_base_java_inputs.py \
  build/sample-singularity-base-source-profile-2026-10-15/stages/n4/attempt-01/t1-n4.nii.gz \
  build/sample-singularity-base-source-profile-2026-10-15/stages/toads-first/attempt-01/dura-inputs/t1-n4_wmfill.nii.gz \
  build/sample-singularity-base-source-profile-2026-10-15/stages/toads-first/attempt-01/dura-inputs/t1-n4_gm.nii.gz \
  build/sample-singularity-base-source-profile-2026-10-15/stages/toads-first/attempt-01/dura-inputs/t1-n4_wmmask.nii.gz \
  build/sample-singularity-base-source-profile-2026-10-15/java-dura/input.bin

PYTHONPATH=.:src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- \
  env DURA_WRAPPER_NX=184 DURA_WRAPPER_NY=256 DURA_WRAPPER_NZ=256 \
  DURA_WRAPPER_ITERATIONS=5 DURA_WRAPPER_MEMBERSHIP=0.3 \
  DURA_WRAPPER_GM_REGULARIZE=0.02 DURA_WRAPPER_JAVA_HEAP=8g \
  scripts/run_java_dura_base_reference.sh \
  build/sample-singularity-base-source-profile-2026-10-15/java-dura/input.bin \
  build/sample-singularity-base-source-profile-2026-10-15/java-dura/output
```

Then run `scripts/compare_dura_java_outputs.py` with the same four NIfTIs, the
Java output directory, `--membership-level 0.3`, and `--iterations 5`. The
machine-readable comparison and intermediate outputs are under the ignored
`build/sample-singularity-base-source-profile-2026-10-15/java-dura/` directory.

## Interpretation

Combined with the exact held-fixed Java matches for both TOADS handoffs in the
parent report, these results establish selected-source parity for the two
TOADS classifications and the intervening Dura stage on this human-derived
Python workflow. They do not resolve differences between the downloaded
historical output and this run: its original effective settings and source
provenance are unavailable, and the upstream N4/brain-mask contract is not
matched.
