# Full pytest report — LesionTOADS output I/O — 2026-10-01

## Result

- Tested source tree, committed unchanged afterward as: `006482f03a28809b0cf01b92b2f35f6f7b9f9111`
- Command: `./.venv/bin/pytest -q`
- Result: **1,028 passed, 0 failed, 0 skipped**
- Expected RAS-orientation warnings: **11**
- Elapsed time: **1,290.93 seconds** (21 minutes, 30.93 seconds)
- Data scope: committed synthetic fixtures and bundled fixed atlas resources only; no subject imaging data.

## Warning categories

The 11 warnings are expected notices for non-RAS synthetic inputs: two image-calculator inputs, one N4 input, two Dura-topology inputs, and six workflow-preflight inputs. No warning indicated a failed test.

## Focused evidence

The related P39 output adapter, existing output dispatch/route tests, plus-input handoff, and image I/O checks passed together: **15 passed in 74.36 seconds**. The new adapter tests additionally verify unsigned-byte label preservation, signed-short subtype labels, float32 channels, class/member 4-D dimensions, RAS affine handling, `.nii.gz` paths, collision refusal, and wrong-grid rejection.

## Scope and limits

This validates conversion and serialization of the already-uncropped selected LesionTOADS output arrays, using the RAS T1 image for physical geometry. It does not establish MIPAV `ModelImage` lifecycle equivalence, full LesionTOADS execution, or complete MaCRUISE workflow equivalence. Synthetic-only evidence is preserved; the supplied human MRI run was not used.
