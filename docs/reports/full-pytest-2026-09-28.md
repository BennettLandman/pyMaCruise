# Full pytest evidence for completed roadmap rows

This report records full-repository test runs for rows marked **Done** in
`plan.md`. Each row identifies its stage-specific evidence and cites the exact
full-suite run and tested implementation commit. Full-suite results provide
regression evidence; they do not replace the row-specific numerical
comparisons, source audits, asset checks, and decisions documented in the
linked evidence.

## Recorded full-suite run

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **212 passed in 71.33 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `0852b56c4dbf94dbd973c31a819110426316652b` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 212 tests collected and all passed |

The tested commit is the implementation and test state used for this run.
Only roadmap and report documentation was updated afterward to record its
results and add P18i; implementation and test files are unchanged.
The per-row entries identify this run wherever it is the cited validation
evidence; later rows can cite a newer recorded run.

## Additional full-suite run for P18j

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **218 passed in 63.57 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `c06d00c654c7642b708049aeeb95c221c35b1dfb` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 218 tests collected and all passed |

This run adds the P18j local TOADS multi-label topology gate. Earlier row
entries retain their original full-suite run and exact tested commit; P18j
cites this newer run.

## Additional full-suite run for P18k

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **223 passed in 48.40 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `2ec8971da2c6e012c0cafb9f0669e680f043868b` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 223 tests collected and all passed |

This run validates the P18k transformed shape-support bounds slice, including
the follow-up that keeps per-voxel support-coordinate recording opt-in so
normal atlas scans retain only bounds and a count. Earlier row entries retain
their respective full-suite runs and exact tested commits; P18k cites this
latest run.

## Additional full-suite run for P18m

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **237 passed in 69.89 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `365ef755c7f0cfb5a927001b02f0d8ae68fe3bd4` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 237 tests collected and all passed |

This run validates the isolated P18m exact-membership registration energy
kernel. It does not validate rigid transform estimation, optimizer convergence,
TOADS evolution, or pipeline equivalence. Earlier row entries retain their
respective full-suite runs and exact tested commits; P18m cites this latest run.

## Additional full-suite run for P18n

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **239 passed in 68.34 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `57c2f42dc767ef901802ab6898656320622c9049` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 239 tests collected and all passed |

This run validates the P18n source call-order excerpt, including the
source-faithful coefficient normalization behavior and the opt-in corrected
calculation. It does not validate coefficient accumulation or registration
quality. Earlier rows retain their respective full-suite evidence; P18n cites
this run.

## Additional full-suite run for P18o

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **240 passed in 68.24 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `41e6e0cad3b9bd17dd0b982fcd5fd0a27a2a69d2` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 240 tests collected and all passed |

This run validates the P18o trilinear prior-derivative kernel on four synthetic
coordinates against a source-formula Java excerpt, including boundary cases.
It does not validate rigid-transform derivatives or registration coefficients.
Earlier rows retain their respective full-suite evidence; P18o cites this run.

## Additional full-suite run for P18p

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **242 passed in 69.63 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `8e7a0c070f325349a5f24b006864bda46c26f10b` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 242 tests collected and all passed |

This run validates the P18p rigid coordinate Jacobian on the recorded
non-singular synthetic case and rejects a singular unit-norm derivative case.
It does not validate prior-gradient composition or registration quality.
Earlier rows retain their respective full-suite evidence; P18p cites this run.

## Additional full-suite run for P18q

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **243 passed in 69.40 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `309f7d66c672fb9c7d841928f62127a93bc3bc4d` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 243 tests collected and all passed |

This run validates P18q's one-voxel exact coefficient contributions against a
source-formula Java excerpt. It does not validate full-volume accumulation,
coefficient normalization, optimizer behavior, or registration quality.
Earlier rows retain their own full-suite evidence; P18q cites this run.

## Additional full-suite run for P18r

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **244 passed in 68.83 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `f44e416f4fab1d36d36583ea336385ed107dac63` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 244 tests collected and all passed |

This run validates P18r's one-step source line search on synthetic accepted
and rejected energy callbacks, including the selected 11-rejection stop. It
does not validate integrated energy/gradient evaluation or atlas alignment.
Earlier rows retain their own full-suite evidence; P18r cites this run.

## Additional full-suite run for P18s

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **245 passed in 68.98 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `5ba9abc87852f77787748a005c4586c1c3bdd36f` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 245 tests collected and all passed |

This run validates P18s full sampled-grid coefficient accumulation on two
small synthetic scans against the Java source-formula excerpt. It does not
integrate transform-dependent atlas sampling or establish large-volume
performance. Earlier rows retain their own full-suite evidence; P18s cites
this run.

## Evidence by Done row

| Plan row | Row-specific validation evidence | Full-suite evidence |
|---|---|---|
| P01 — Establish this roadmap | `README.md` links the roadmap; `plan.md` contains ordered criteria and status/date fields. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P02 — Select canonical base workflow | `docs/source-baseline.md` records the V3 source anchor and version clarification. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P03 — Select supported MaCRUISEplus contract | `docs/plus-baseline.md` documents the selected profile and deferred variants. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P04 — Build source-to-stage traceability map | `docs/source-map.md` maps the selected source stages and gaps. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P05 — Inventory external code and assets | `docs/assets-and-dependencies.md` records source trees, revisions, dependencies, and missing assets. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P05a — Bundle supplied workflow atlases | `tests/test_atlas_assets.py` validates the bundled atlas inventory and image assets. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P05b — Preserve license and attribution information | `license_info.md` and `licenses/upstream/` contain the preserved notices and handoff overview. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P05c — Normalize bundled image assets to RAS `.nii.gz` | `tests/test_atlas_assets.py` checks normalized atlas payloads and descriptors; audit details are in `data/atlases/README.md` and `docs/assets-and-dependencies.md`. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P06 — Freeze image and coordinate conventions | `tests/test_image_io.py`, `tests/test_surfaces.py`, `tests/test_spatial_transfer.py`, plus `docs/image-coordinate-conventions.md`. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P07 — Define workflow input/output contracts | `tests/test_contracts.py` and `docs/io-contract.md`. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P08 — Define numerical acceptance policy | `tests/test_validation.py` and `docs/validation-tolerances.md`. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P09 — Create deterministic synthetic phantom generator | `tests/test_synthetic.py` validates seeded fixture generation and metadata. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P10 — Create synthetic fixture corpus | `tests/data/synthetic/` inventory and `tests/test_synthetic.py`; atlas fixtures additionally use `tests/test_atlas_assets.py`. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P11 — Establish test and comparison harness | `tests/test_validation.py` validates numerical and geometric comparison helpers. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P12 — Create packaging and dependency profiles | `tests/test_cli.py` checks installed command behavior; package and dependency declarations are in `pyproject.toml`. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P13 — Implement geometry-aware image I/O | `tests/test_image_io.py` validates RAS normalization, round-trip behavior, and failure cases. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P14 — Implement configuration, logging, and run manifest | `tests/test_workflow.py` and `src/pymacruise/workflow.py` cover manifest/configuration behavior. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P15 — Validate existing generic volume operator | `tests/test_volume_ops.py` validates array operations, geometry, dtype, and NIfTI behavior. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P16 — Port raw input preparation | `tests/test_raw_preparation.py` validates synthetic RAS preparation behavior. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P17 — Port cortical separation and label preparation | `tests/test_cortical_separation.py` validates source label sets and output masks. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18a — Port TOADS robust intensity normalization | `tests/test_toads_preprocessing.py`, `tests/java_reference/ToadsIntensityReferenceHarness.java`, and [docs/toads-port-plan.md](../toads-port-plan.md) record exact synthetic Java matching. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18b — Implement TOADS initial centroids and exact memberships | `tests/test_toads_initialization.py`, `tests/java_reference/ToadsInitializationReferenceHarness.java`, and `tests/data/synthetic/toads_initialization/` compare all centroid and membership values against the Java harness, including tied-centroid weighting and the zero-distance cap. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18c — Parse selected TOADS atlas intensity priors | `tests/test_toads_atlas_priors.py` and its Java descriptor harness verify the selected base atlas's T1 MPRAGE class priors, class count/order, aliases, channel order, and malformed/missing profile handling. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18d — Port original-JIST dura signed distance field | `tests/test_toads_dura.py` compares the adapter's float32 outputs with all four original-JAR fields in `tests/data/synthetic/p18_distancefield_jar_reference/`; 648/648 voxels match exactly, and invalid input cases reject. This is limited to signed-distance construction. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18e — Port TOADS crop and topology geometry | `tests/test_toads_crop.py` compares bounds, padded image arrays, transform changes, uncropped arrays, and categorical topology with the documented formula-transcription fixture. Bounds/labels/float32 image values match exactly; translations are within `2e-7`. This is not an upstream class execution. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18f — Bundle selected TOADS 18/6 critical-point LUT | `tests/test_toads_critical_lut.py` validates the runtime packaged resource, payload length/domain, read-only array, malformed inputs, and Java x/y/z key packing. `tests/test_critical_lut_source.py` checks source SHA, gzip size, representative keys, and the known distinction from P25's 6/18 predicate. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18g — Port Dura dynamic-pressure point update | `tests/test_toads_dura.py` compares an interior positive-pressure point and a mirrored-boundary negative-pressure point against the selected VandyCRUISE `DynamicPressureTGDM.updatePoint` and `DynamicPressureCortexProfile` Java harness within `5e-10`; invalid shapes and coordinates reject. This does not validate the GenericTGDM evolution loop, topology veto, scheduling, or reinitialization. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18h — Use selected TOADS 18/6 LUT predicate | `tests/test_toads_critical_lut.py` compares key packing and lookup result for all 68 committed neighborhoods against the selected Java `CriticalPointLUT` class (four edge cases plus 64 seeded patterns). This validates one-neighborhood regularity only, not multi-label homeomorphism checks. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18i — Port fixed-transform TOADS shape sampling | `tests/test_toads_shape_sampling.py` compares eight fixed-transform prior samples, mapped coordinates, boundary defaults, and slope/clipping against the selected Java `ParametricTransform`, `ImageFunctions`, `Numerics`, and shape postprocessing harness. Float32 outputs and mapped coordinates match exactly. This does not run atlas registration or optimized bounds. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P18j — Port TOADS 18/6 multi-label homeomorphic-change gate | `tests/test_toads_topology.py` compares 130 synthetic cases against the selected Java method-body excerpt and actual `CriticalPointLUT`; fixture scope is documented under `tests/data/synthetic/toads_homeomorphic_change/`. | `218 passed`; commit `c06d00c654c7642b708049aeeb95c221c35b1dfb` |
| P18k — Port transformed TOADS shape support bounds | `tests/test_toads_shape_bounds.py` compares all positive-support coordinates and bounds in five cases against the Java transform/interpolation methods plus selected atlas source scan excerpt; also verifies max-exclusive and empty-support semantics. The production call retains only bounds and support count. | `223 passed`; commit `2ec8971da2c6e012c0cafb9f0669e680f043868b` |
| P18m — Port exact-membership rigid-registration energy kernel | `tests/test_toads_registration.py` compares two synthetic objective accumulations against the Java source-formula harness and checks source empty-norm NaN behavior and invalid inputs. This is kernel evidence only, not registration equivalence. | `237 passed`; commit `365ef755c7f0cfb5a927001b02f0d8ae68fe3bd4` |
| P18n — Preserve TOADS registration coefficient normalization call order | `tests/test_toads_registration.py` matches the fixed source-formula harness in default source-faithful mode and separately checks that the apparent corrected gradient is opt-in. | `239 passed`; commit `57c2f42dc767ef901802ab6898656320622c9049` |
| P18o — Port TOADS trilinear prior derivatives | `tests/test_toads_registration.py` compares x/y/z prior derivatives on four valid and out-of-bounds synthetic coordinates against the source-formula Java excerpt. | `240 passed`; commit `41e6e0cad3b9bd17dd0b982fcd5fd0a27a2a69d2` |
| P18p — Port TOADS rigid transform coordinate Jacobian | `tests/test_toads_registration.py` compares the 3×6 rotation/translation Jacobian against the Java formula excerpt and checks singular-rotation rejection. | `242 passed`; commit `8e7a0c070f325349a5f24b006864bda46c26f10b` |
| P18q — Port per-voxel TOADS registration coefficient accumulation | `tests/test_toads_registration.py` compares cost, norm, gradient, diagonal Hessian, norm-gradient, and norm-Hessian values for a synthetic voxel against a Java method-formula excerpt. | `243 passed`; commit `309f7d66c672fb9c7d841928f62127a93bc3bc4d` |
| P18r — Port one TOADS registration gradient-descent line-search step | `tests/test_toads_registration.py` compares accepted/rejected parameter states, lambda, relative change, and retry count against the Java source-formula excerpt. | `244 passed`; commit `f44e416f4fab1d36d36583ea336385ed107dac63` |
| P18s — Port full sampled-grid TOADS registration coefficient accumulation | `tests/test_toads_registration.py` compares full-grid and offset/subsample source-order float32 sums to a Java formula excerpt. | `245 passed`; commit `5ba9abc87852f77787748a005c4586c1c3bdd36f` |
| P19 — Integrate ANTsPy N4 bias correction | `tests/test_n4.py` checks synthetic bias reduction, geometry, and repeatability. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P20 — Port crop and inverse-crop transforms | `tests/test_cropping.py` validates synthetic crop bounds and inverse geometry. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P21 — Port TOADS outputs back to source space | `tests/test_spatial_transfer.py` and the Java/Python comparison artifacts cited in `plan.md`. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P22 — Port multi-atlas boundary/CSF generation | `tests/test_boundary.py` validates boundaries, exclusions, and empty cases. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P23 — Port boundary cleanup | `tests/test_boundary.py` validates cleanup behavior on designed synthetic arrays. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P24 — Port multi-atlas and TOADS tissue combination | `tests/test_boundary.py` validates priority, erosion, and grid behavior. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P25 — Port topology correction | `tests/test_topology.py` validates local connectivity, components, cavities, and repeatability. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P26 — Port ACE tissue enhancement | `tests/test_ace.py`, `tests/test_ace_reference.py`, and `docs/ace-port.md` provide Java-stage and numerical evidence. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P27 — Port enhanced skeleton generation | `tests/test_skeleton_reference.py` and `docs/skeleton-port.md` record source comparisons and surface metrics. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P28 — Port Gradient Vector Flow | `tests/test_gvf.py` and `docs/gvf-port.md` record the compiled-Java synthetic comparison. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P29 — Port CRUISE/TGDM level-set evolution | `tests/test_tgdm.py`, `docs/tgdm-reference-comparison.md`, and the committed Java-reference JSON report. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P30 — Port level-set to surface extraction | `tests/test_surfaces.py` and `docs/surface-extraction.md` document the analytic mesh checks and Lewiner substitution. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P31 — Port scientific surface and thickness outputs | `tests/test_thickness.py` and `docs/thickness-port.md` provide known-gap and geometry results. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P32 — Port segmentation refinement | `tests/test_segmentation_refine.py` and `docs/segmentation-refinement.md` record synthetic threshold and rule evidence. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P33 — Port MNI registration and transfer | `tests/test_registration.py`, `tests/test_spatial_transfer.py`, and `docs/registration-port.md` record synthetic transform evidence. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P38 — Implement MaCRUISEplus FLAIR registration | `tests/test_registration.py` and `docs/registration-port.md` record the known-shift synthetic test and configuration checks. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P39a — Port LesionTOADS intensity and centroid initialization | `tests/test_lesion_toads_preprocessing.py`, its test-only Java harness, and the committed two-channel RAS `.nii.gz` phantom compare robust bounds, all normalized values, and lesion/black-hole centroids; zero-span rejection is documented as a deliberate invalid-input deviation. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P39b — Parse LesionTOADS atlas priors | `tests/test_lesion_atlas_priors.py` and the image-free descriptor fixture verify structure order, selected T1 MPRAGE/FLAIR rows, lesion/black-hole priors, source aliases, and malformed/incomplete descriptor handling. The bundled atlas values match the selected source descriptor. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |
| P39c — Port LesionTOADS capped distance support | `tests/test_lesion_toads_distance.py` compares all synthetic mask fields against a standalone transcription of the selected VandyCRUISE boolean-mask helper; per-voxel exact equality includes the source's cap-1 behavior. This is not lesion segmentation equivalence. | `212 passed`; commit `0852b56c4dbf94dbd973c31a819110426316652b` |

## Future rows

When a roadmap row is completed, run the complete suite and add its row-specific
evidence, full-suite pass count, exact tested commit, and run date here in the
same change that marks the row Done. If implementation changes after a recorded
run, rerun the complete suite before adding or refreshing completion evidence.
