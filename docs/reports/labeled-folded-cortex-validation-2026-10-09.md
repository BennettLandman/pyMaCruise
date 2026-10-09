# Labeled folded-cortex TGDM validation

**Implementation commit:** `eb04598` (`Add labeled folded cortex TGDM validation case`); saved-stage datatype correction is included in the follow-up report commit.  
**Evidence categories:** synthetic integration passed; full-profile selected Java surface gates passed; no human MRI used.

## Case and runnable checks

The committed `tests/data/synthetic/labeled_folded_cortex/` case is a 65³,
1 mm isotropic RAS image with MaSeg-like labels 0 (background), 5 (sulcal
CSF), 44 (WM), and 100 (cortical GM). The folded shell has an analytic nominal
5 mm ribbon. Gray-matter support is 255 inside the ribbon; the optional ACE
skeleton is zero to isolate the nested TGDM surfaces. Inputs, generator,
assumptions, and Java references are documented beside the fixture.

Run the complete selected profile under an 8 GB sampled-RSS watchdog:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/run_labeled_folded_cortex.py \
  --output-dir build/labeled-folded-cortex-source-selected
```

The run used full selected iterations: smoothing 6, initial/inner 4,
central 7, central shrink 1, and outer 7 (the selected wrapper maps the
central iteration option to outer). It used WM iso 0.6 and zero GVF. It
completed successfully with peak sampled RSS **0.336 GB**. The binary mask is
converted to the inside-negative/outside-positive sign field expected by
`GenericTGDM`'s bound logic; passing raw 0/255 mask values blocked evolution.
The adapter correction is covered by two mask tests, and the full focused
TGDM module passed **26 tests**.

Synthetic ribbon checks passed:

- All input arrays are finite and share the RAS grid.
- Inner, central, and outer meshes are closed, manifold, connected, and have
  Euler characteristic 2; voxel regions are nested with zero violations.
- WM is inside inner at **100%**, GM is between inner and outer at **99.983%**,
  and CSF is outside outer at **100%**. Four GM voxels fell outside outer.
- Symmetric nearest-vertex inner/outer separation is finite and positive;
  median **4.757 mm**, 5th percentile **4.341 mm**, 95th percentile **5.091
  mm**. This is a sampled diagnostic, not exact point-to-triangle thickness.

The detailed machine-readable report is
[labeled-folded-cortex-ribbon-quality-2026-10-09.json](labeled-folded-cortex-ribbon-quality-2026-10-09.json).
It is reproducible from the command above.

## Selected Java comparison

Java references were generated with the supplied `MaGenericTGDM`/`GenericTGDM`
and `DistanceField` validation sources on the same NIfTIs, with WM iso 0.6,
the explicit signed brain bound, and full settings. The Java run used **0.597
GB peak sampled RSS** under an 8 GB watchdog. Reproduce it with
`scripts/run_labeled_folded_cortex_java.sh`; compare with
`scripts/compare_labeled_folded_cortex_java.py`.
The compiled reference inputs are pinned by SHA-256: `GenericTGDM.java`
`70f5b215…`, `MaGenericTGDM.java` `66124d41…`, `DistanceFieldyh.java`
`a0df8859…`, the omitted-class `DistanceField` bridge `91c0e8b9…`, and
`TGDMReferenceHarness.java` `a8114822…`.

All three final regions passed all agreed gates (per-region Dice ≥ 0.97,
symmetric mean distance ≤ 0.01 mm, Hausdorff ≤ 1 mm):

| Surface | Dice | Mean distance (mm) | Hausdorff (mm) | Result |
| --- | ---: | ---: | ---: | --- |
| Central | 1.000000 | 0.000045 | 0.017552 | Pass |
| Inner | 1.000000 | 0.00000014 | 0.000004 | Pass |
| Outer | 1.000000 | 0.000773 | 0.070323 | Pass |

The selected Java/Python level-set MAE is 2.8e-8 for the signed seed, 1.31e-4
after initial smoothing, 8.39e-6 for initial-inner, 4.73e-4 for central,
1.25e-4 for central-shrink, 2.43e-7 for final inner, and 4.97e-4 for outer.
The exact Java/Python stage metrics are in
[labeled-folded-cortex-java-comparison-2026-10-09.json](labeled-folded-cortex-java-comparison-2026-10-09.json).

The initial comparison attempt was discarded: it copied the uint8 input NIfTI
header onto float level sets, quantizing saved Python files. The runner now
sets and verifies each saved NIfTI datatype (float32 for level sets, uint8 for
tissue labels); the comparison and JSON report were regenerated from those
verified float32 outputs. No threshold was changed.

## Related TOADS evidence

The corrected second TOADS pass's `DuraRemoval_input` handoff was independently
replayed with the selected Java wrapper from the exact saved Python input. The
five products (classification, CSF, GM, WM fill, and WM mask) match exactly
after the wrapper's unconditional ×255 scale; all five products also matched
exactly in the first pass. The second-pass case is 172×216×164 RAS, naturally
stopped at 43 iterations in the integrated workflow; this source replay used
its saved handoff arrays. Java peak sampled RSS was **3.722 GB** under a 16 GB
watchdog. See
[toads-integrated-second-dura-java-comparison-2026-10-09.json](toads-integrated-second-dura-java-comparison-2026-10-09.json).

This confirms the Dura-to-TOADS handoff products on those held-fixed inputs; it
does not establish full MaCRUISE workflow equivalence or real-data validation.
