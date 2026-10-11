# In-vivo TOADS Java source localization — 2026-10-10

## Result

The supplied human T1 has a direct selected-Java TOADS comparison. It localizes the discrepancy after initialization and before natural convergence:

- With the same masked, N4-corrected first-pass T1, Java and Python zero-iteration classifications are **byte-for-byte identical** on the `(208, 309, 305)` RAS grid (19,602,960 voxels). Both raw maps have SHA-256 `d38fc4f7168b3fee81b762bf22b3f851ae8c430b8d88c72186a1f7ae19f9c9af`.
- The four-level registration and post-refresh transforms also match Java at the first pass, including recropped translation `(-4.2854285, -1.3406382, -7.2076693)` mm. Atlas topology and all 11 shape maps were checked against the Java raw assets and were value-exact.
- After the first TOADS iteration, Java/Python differ at **3,234 voxels** (0.0165%), with minimum label Dice **0.978335** (label 6); all labels still pass the 0.97 threshold.
- After two TOADS iterations on that same first-pass input, Java/Python differ at **7,057 voxels** (0.0360%). Every label Dice is at least **0.980025** (minimum: label 6); the maps do not match exactly, but all pass the 0.97 label threshold at this short checkpoint.
- At natural convergence, the full source replay still fails the per-label threshold. Java/Python stopped at iterations **23/22** on first pass and **24/22** on second pass. The first-pass Dice minima are label 6 **0.954809** and label 16 **0.966437**; second-pass minima are label 6 **0.966144** and label 16 **0.966268**. Full voxel differences are 10,809 and 13,762, respectively.

This rules out the atlas conversion, initial Java/Python class assignment, and first-pass rigid registration as the source of the observed gap. It puts the first observed difference inside the first evolving cycle; the disagreement grows across later cycles until two required labels fail. It does not yet isolate which membership, thinning, growth, relation, or registration update introduces it.

## First-iteration label results

| Label | Java voxels | Python voxels | Dice |
|---:|---:|---:|---:|
| 1 | 18,254,941 | 18,254,827 | 0.999991 |
| 5 | 232,511 | 232,631 | 0.996182 |
| 6 | 3,972 | 3,921 | 0.978335 |
| 14 | 101,200 | 101,167 | 0.998809 |
| 15 | 609,726 | 609,617 | 0.998338 |
| 16 | 3,393 | 3,343 | 0.972981 |
| 17 | 15,756 | 15,733 | 0.997047 |
| 18 | 7,937 | 7,958 | 0.997169 |
| 22 | 20,356 | 20,381 | 0.997864 |
| 24 | 28,607 | 28,600 | 0.996102 |
| 25 | 324,561 | 324,782 | 0.998034 |

The Java/Python first-iteration raw output hashes are `713ed9427c1d3572ed4101af07e04e8fad18ad704721356b4b7e89999fae4651` and `040d3228134fcd76a517b206d5a0af10ee79624bc7c85cbf471ab1ce34b8e0b3`, respectively.

## Two-iteration label results

| Label | Java voxels | Python voxels | Dice |
|---:|---:|---:|---:|
| 1 | 18,266,324 | 18,266,560 | 0.999950 |
| 5 | 184,291 | 183,795 | 0.986297 |
| 6 | 4,368 | 4,393 | 0.980025 |
| 14 | 109,870 | 109,933 | 0.997493 |
| 15 | 638,194 | 638,462 | 0.996980 |
| 16 | 3,893 | 3,894 | 0.986516 |
| 17 | 16,478 | 16,477 | 0.994082 |
| 18 | 8,198 | 8,214 | 0.997075 |
| 22 | 20,152 | 20,204 | 0.992467 |
| 24 | 30,742 | 30,695 | 0.995588 |
| 25 | 320,450 | 320,333 | 0.997308 |

At natural convergence, the full first- and second-pass label tables are retained in ignored `build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/first-comparison-masked.json` and `second-comparison-masked.json`.

## Method and reproducibility

The Java reference uses the test-only `ToadsWorkflowSourceHarness` and the original VandyCRUISE text descriptor plus raw atlas assets. Its two-iteration invocation is:

```sh
CLASSPATH=build/java-reference/toads-selected-wrapper-harness:$(cat build/java-reference/toads-selected-wrapper-harness/classpath.txt)
java -Xmx12g -cp "$CLASSPATH" ToadsWorkflowSourceHarness \
  build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/first-input-masked.raw \
  /Users/landmaba/dev/mipav/toads-cruise-scm-2026-09-26/VandyCRUISE/src/Atlas/cruise-atlas-10obj-toads2008.txt \
  208 309 305 1 1 1 2 \
  build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/first-classification-java-two-iterations.raw
```

The matching Python calls used `register_toads_from_atlas_assets` with four registration levels, 50 registration iterations, and five post-refresh registrations, then `run_registered_toads_loop` with the source profile and a two-iteration ceiling. Classification arrays were compared in Java's x-fastest raw order. The exact-equal zero-iteration check used the same pipeline and settings with zero evolution iterations.

Relevant input atlas conversion audit confirms byte/value equality for topology and all 11 shape maps. The harness is committed in PyMaCRUISE commit `491869ba71214cd7fb931f86c8df449660b67ad4`; its SHA-256 is `397cb71fc3a8a18186fc132321bae54b7aa3691a9fd5331725f1b8287885c2c2`. Selected Java source hashes are: `AlgorithmMultipleToads.java` `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`, `MultipleToadSegmentation.java` `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`, and the original descriptor `d43efa5a002b0c1a838449ba7ca1d3bf35fab5e3b575474a71aadfe475dcc257`. Human images, raw arrays, intermediate products, and logs remain under ignored `build/` and were not committed.

This is **source comparison on one in-vivo case**, not whole-MaCRUISE equivalence or clinical validation. It does not compare final surfaces with historical Java/MATLAB products, and it does not pass the agreed natural-convergence Dice requirement for labels 6 and 16. Next: capture and compare selected intermediate states at the first evolving iteration, then repeat the converged Java/Python comparison after any correction. Source-line delta remains **0** because this investigation only connects already-counted methods.
