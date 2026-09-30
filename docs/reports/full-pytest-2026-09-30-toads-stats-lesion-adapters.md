# TOADS statistics and LesionTOADS output adapters — 2026-09-30

## Full-suite evidence

- Command: `.venv/bin/python -m pytest -q`
- Tested implementation and fixture-provenance commit: `0c9af567c7956bb6886bedb570b4fa50bdede6c2`
- Result: **681 passed, 2 expected RAS reorientation warnings, 158.56 seconds**.
- Both warnings are emitted by `test_dura_first_wm_topology_call_uses_source_adapter_settings_and_ras`; that test deliberately supplies LPS-oriented topology and paint-mask inputs and checks their conversion to RAS.
- Focused command: `.venv/bin/python -m pytest -q tests/test_lesion_gm_ventricle_membership.py tests/test_lesion_cruise_inputs.py tests/test_lesion_dura_outputs.py tests/test_toads_object_centroids.py tests/test_toads_object_variance.py tests/test_toads_object_image_variance.py` — **17 passed in 2.32 seconds**.

This report supplies the full-suite evidence for P18cq–P18cs and P39zzzd–P39zzzzf. All new image fixtures are synthetic NIfTI files with RAS axis codes; the Java references are method-excerpt harnesses, not integrated workflow references.

## P18cq — ordinary-OBJECT centroid update

The selected VandyCRUISE `MultipleToadSegmentation.computeExactCentroids`
method (source SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`)
updates ordinary OBJECT centroids with membership-squared accumulation over
the one-voxel interior. The adapter also matches the active damped mode and
wrapper smoothness `0.1`; cases cover optional field multiplication, PSF
ridge boundary Laplacian, and zero-support retention. Focused module: 2 passed.
This does not implement centroid dispatch for other object types or a full
TOADS iteration.

## P18cr — ordinary-OBJECT assigned-class variance

The selected `computeExactVariances` ordinary OBJECT branch scans the full
volume, applies optional field correction, weights squared residuals by
`membership² + 1e-30f`, and returns zero for classes with no assigned voxels.
The adapter reuses the established float32 variance kernel rather than
duplicating it. Java comparisons cover boundary-only assigned support,
empty support, and zero memberships. Focused module: 4 passed. Lesion-specific
variance splitting and approximate variance are outside this slice.

## P18cs — ordinary-OBJECT image-wide variance

The selected `computeExactImageVariance` ordinary OBJECT branch is exposed
through a class-type-guarded adapter to the existing image-wide variance
kernel. The Java fixture checks raw accumulators, the `1e10f` zero-numerator
fallback, normalized inverse-variance modality weights, standard deviations,
field correction, and full-grid boundary memberships. Focused module: 4
passed. This does not implement the exact membership estimator or special
LESION/OUTLIER mixtures.

## P39zzzd — ordinary GM ventricle-distance attenuation

The selected VandyCRUISE `LesionToadSegmentation.computeMemberships` branch
attenuates only exact ordinary class names `Cerebrum-GM` and `GM`, only when
lesions exist and ventricle distance is valid with nonzero maximum. DIRECT
and PRIOR expression order, zero-distance behavior, fallthroughs, and
positive-denominator normalization match Java float bits. Focused module: 2
passed. Upstream class terms and complete LesionTOADS iteration are not
included.

## P39zzzze — LesionTOADS CRUISE input channels

The selected VandyCRUISE `LesionToadSegmentation.generateCruiseInputs`
method (source SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`)
produces CSF, GM, WM, and hard-WM channels. Synthetic comparisons preserve
exact aliases, last-match CSF/GM overwrite, WM accumulation exclusions, and
hard-WM selection. Focused module: 2 passed. This is an input adapter, not
CRUISE or lesion-segmentation equivalence.

## P39zzzzf — LesionTOADS dura-removal output channels

The selected VandyCRUISE `LesionToadSegmentation.generateDuraRemovalOutputs`
adapter is separate from the base TOADS Dura output routine. Tests cover
CSF overwrite, summed cortical/cerebellar GM, WM fill, hard-WM behavior, and
the historical `CerebellarGM)` spelling predicate. The four outputs match
the pinned excerpt on synthetic RAS NIfTIs. Focused module: 3 passed. This
does not implement the subsequent dura-removal segmentation stages.
