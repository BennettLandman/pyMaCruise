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

## Additional full-suite run for P18t

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **248 passed in 78.71 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `e45528ed129df5628f505f288a7281d73e4a3849` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 248 tests collected and all passed |

This run validates P18t's dense transform-dependent registration-field
composition on small synthetic grids. Existing focused reference tests
independently compare the scalar transform/interpolator, trilinear derivative,
and rigid coordinate Jacobian kernels against copied Java APIs or source-formula
references. It does not execute integrated coefficient optimization; the dense
field helper is not suitable for production-size volumes, where the next
registration driver must stream voxel contributions. Earlier rows retain their
own full-suite evidence; P18t cites this run.

## Additional full-suite run for P18u

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **251 passed in 78.79 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `bdc233277c2227049d6fe969cefbe7c929e5f983` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 251 tests collected and all passed |

This run validates the small-grid TOADS rigid optimizer integration and its
deterministic repeat/stop behavior. It is not a full Java registration
comparison or alignment-quality assessment; the dense implementation is not
ready for clinical-size grids. Earlier rows retain their own full-suite
evidence; P18u cites this run.

## Additional full-suite run for P18v

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **252 passed in 70.82 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `0fe1161d6d5cf8f4af61b13346801ce9236d49dc` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 252 tests collected and all passed |

This run validates the streamed energy and coefficient evaluator against the
dense small-grid path, optimizer repeatability/stops through streaming, and
the synthetic scratch-memory/time benchmark. Results and limitations are in
[the P18v report](toads-registration-streaming-2026-09-28.md). Earlier rows
retain their own full-suite evidence; P18v cites this run.

## Additional full-suite run for P18w

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **253 passed in 75.48 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `dad38b855a480638eb717c0d90421b5c8aa9de38` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 253 tests collected and all passed |

This run validates the vectorized source-order registration kernels against
both scalar streaming and dense synthetic references and reruns the optimizer
checks through the vectorized path. The performance figures and limits are in
[the P18w benchmark report](toads-registration-vectorization-2026-09-28.md).
Earlier rows retain their own full-suite evidence; P18w cites this run.

## Additional full-suite run for P18x

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **259 passed in 50.60 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `d855faf8e49aab850eb4c47e71ec7fd14702651f` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 259 tests collected and all passed |

This run validates the source-order pyramid block means and coarse-to-fine
scale schedule on synthetic arrays. Exact reference values and remaining
scaled-coordinate work are recorded in [the P18x report](toads-registration-pyramid-2026-09-28.md).
Earlier rows retain their own full-suite evidence; P18x cites this run.

## Additional full-suite run for P18y

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **263 passed in 57.24 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `0854326892f347b95b58d7d582f116674293afa3` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 263 tests collected and all passed |

P18y adds exact synthetic comparisons for coarse-scale rigid coordinates and
all 18 entries of the rigid coordinate Jacobian at scales 1, 2, and 4. Earlier
rows retain their respective full-suite evidence; P18y cites this run.

## Additional full-suite run for P18z

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **266 passed in 48.40 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `2b42a9567e6cbdfc4a59da316f87ba28ba838df2` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 266 tests collected and all passed |

P18z adds coarse-to-fine rigid optimizer integration at source scales 4 and 2,
including committed synthetic RAS NIfTI inputs, deterministic transform
handoff, the scale-1 no-op contract, and scaled-kernel comparisons. Earlier
rows retain their respective full-suite evidence; P18z cites this run.

## Additional full-suite run for P18aa

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **267 passed in 55.55 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `3eb7215847917d8d6dda0fa97bbc295b3107f436` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 267 tests collected and all passed |

P18aa adds Java-formula reference comparisons for source-scaled registration
energy, norm, gradient, diagonal Hessian, norm-gradient, norm-Hessian, and
positive/zero prior counts at scales 2 and 4. Earlier rows retain their
respective full-suite evidence; P18aa cites this run.

## Additional full-suite run for P18ab

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **268 passed in 51.06 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `ecbebfe5ed27d5c22e85271471dc4e2f22da3e84` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 268 tests collected and all passed |

P18ab adds a larger folded cortex-like RAS NIfTI shape prior and compares its
positive-support count and bounds against the selected Java transform and
interpolation APIs at three rigid transforms. Earlier rows retain their
respective full-suite evidence; P18ab cites this run.

## Additional full-suite run for P18ac

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **268 passed in 50.27 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `32e6d804e56b312da04d02ec658abac57f36493a` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 268 tests collected and all passed |

P18ac adds source-derived geometry behavior to atlas preflight: descriptor
resolution drives the selected Java algorithm, while conflicting NIfTI
metadata remains reported. Earlier rows retain their respective full-suite
evidence; P18ac cites this run.

## Additional full-suite run for P18ad

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **270 passed in 52.01 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `b0dd1fe4aafd9f78b66b137fa0ed04c5630c65fb` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 270 tests collected and all passed |

P18ad compares direct NumPy power and source-index-quantized NumPy power
against 1,042 cases from the actual selected Java `PowerTable` class. The
small normalized likelihood probe shows direct power changes three class
winners, while source-index rounding preserves them. This does not test the
complete TOADS membership update; P18ad cites this run.

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
| P18l — Map TOADS rigid atlas-registration objective and library substitution boundary | `docs/toads-atlas-registration.md` records the selected objective, wrapper call path, source hash, ANTsPy substitution decision, and measured SciPy interpolation difference; later P18m–P18s rows add isolated source-formula checks. | `245 passed`; commit `5ba9abc87852f77787748a005c4586c1c3bdd36f` |
| P18m — Port exact-membership rigid-registration energy kernel | `tests/test_toads_registration.py` compares two synthetic objective accumulations against the Java source-formula harness and checks source empty-norm NaN behavior and invalid inputs. This is kernel evidence only, not registration equivalence. | `237 passed`; commit `365ef755c7f0cfb5a927001b02f0d8ae68fe3bd4` |
| P18n — Preserve TOADS registration coefficient normalization call order | `tests/test_toads_registration.py` matches the fixed source-formula harness in default source-faithful mode and separately checks that the apparent corrected gradient is opt-in. | `239 passed`; commit `57c2f42dc767ef901802ab6898656320622c9049` |
| P18o — Port TOADS trilinear prior derivatives | `tests/test_toads_registration.py` compares x/y/z prior derivatives on four valid and out-of-bounds synthetic coordinates against the source-formula Java excerpt. | `240 passed`; commit `41e6e0cad3b9bd17dd0b982fcd5fd0a27a2a69d2` |
| P18p — Port TOADS rigid transform coordinate Jacobian | `tests/test_toads_registration.py` compares the 3×6 rotation/translation Jacobian against the Java formula excerpt and checks singular-rotation rejection. | `242 passed`; commit `8e7a0c070f325349a5f24b006864bda46c26f10b` |
| P18q — Port per-voxel TOADS registration coefficient accumulation | `tests/test_toads_registration.py` compares cost, norm, gradient, diagonal Hessian, norm-gradient, and norm-Hessian values for a synthetic voxel against a Java method-formula excerpt. | `243 passed`; commit `309f7d66c672fb9c7d841928f62127a93bc3bc4d` |
| P18r — Port one TOADS registration gradient-descent line-search step | `tests/test_toads_registration.py` compares accepted/rejected parameter states, lambda, relative change, and retry count against the Java source-formula excerpt. | `244 passed`; commit `f44e416f4fab1d36d36583ea336385ed107dac63` |
| P18s — Port full sampled-grid TOADS registration coefficient accumulation | `tests/test_toads_registration.py` compares full-grid and offset/subsample source-order float32 sums to a Java formula excerpt. | `245 passed`; commit `5ba9abc87852f77787748a005c4586c1c3bdd36f` |
| P18t — Build transform-dependent atlas fields for TOADS registration | `tests/test_toads_registration_sampling.py` composes class-bound sampling and prior derivatives on synthetic grids; companion kernel tests validate transform/interpolation and Jacobian references. Dense field storage is explicitly limited to small grids. | `248 passed`; commit `e45528ed129df5628f505f288a7281d73e4a3849` |
| P18u — Integrate TOADS sampled-grid rigid registration optimizer | `tests/test_toads_registration_sampling.py` validates deterministic optimizer history, accepted trials, outer stop gates, and no-positive-prior rejection on a synthetic array fixture. Component kernels retain their separate source-reference checks. | `251 passed`; commit `bdc233277c2227049d6fe969cefbe7c929e5f983` |
| P18v — Stream TOADS rigid registration samples and coefficients | `tests/test_toads_registration_sampling.py` asserts exact float32 equality between dense and streamed energy/coefficients on full and offset/subsample synthetic scans; reproducible time/scratch probe: `scripts/benchmark_toads_registration_streaming.py`. | `252 passed`; commit `0fe1161d6d5cf8f4af61b13346801ce9236d49dc` |
| P18w — Vectorize source-ordered TOADS sampling and accumulation | `tests/test_toads_registration_sampling.py` compares vectorized, scalar, and dense energies/coefficients exactly across two anisotropic grids and full/subsampled scans; reproducible 6³/10³/16³ performance and scratch-memory results: `docs/reports/toads-registration-vectorization-2026-09-28.md`. | `253 passed`; commit `dad38b855a480638eb717c0d90421b5c8aa9de38` |
| P18x — Port source TOADS pyramid block downsampling and level schedule | `tests/test_toads_registration_sampling.py` compares factor-2/4 results with the committed standalone Java harness and checks source scale schedule/dimensions for levels 1–5. | `259 passed`; commit `d855faf8e49aab850eb4c47e71ec7fd14702651f` |
| P18y — Port coarse-scale TOADS rigid coordinates and Jacobians | `tests/test_toads_shape_sampling.py` compares mapped coordinates and the complete rigid Jacobian with a committed standalone Java formula fixture at scales 1/2/4, odd image dimensions, anisotropic resolution, and nonzero parameters; scale and center validation are tested. | `263 passed`; commit `0854326892f347b95b58d7d582f116674293afa3` |
| P18z — Integrate coarse-to-fine rigid atlas registration | Committed synthetic inputs `tests/data/synthetic/toads_registration_pyramid_optimizer/*.nii.gz`; `tests/test_toads_registration_sampling.py` verifies RAS geometry and exact seeded values, scales 4→2 updates, deterministic history/transform handoff, level-1 no-op, and scale-2 mapped-kernel comparisons; coefficient maximum absolute difference is `1.49e-8`. | `266 passed`; commit `2b42a9567e6cbdfc4a59da316f87ba28ba838df2` |
| P18aa — Validate coarse-scale registration objective and coefficients against Java formulas | Committed 4D RAS NIfTI fixtures and `ToadsScaledRegistrationKernelsReferenceHarness.java`; `tests/test_toads_registration_sampling.py` compares energy/coefficients/counts at scales 2 and 4 and confirms scalar positive/zero support counts. Maximum absolute difference is guarded at `1.5e-7`. | `267 passed`; commit `3eb7215847917d8d6dda0fa97bbc295b3107f436` |
| P18ab — Validate transformed support-bound refresh on folded cortex-like priors | Committed folded synthetic RAS `.nii.gz` prior; `tests/test_toads_shape_bounds.py` matches support count and all bounds with the Java transform/interpolation harness for identity, translated, and rotated-plus-translated cases; outputs have distinct bounds. | `268 passed`; commit `ecbebfe5ed27d5c22e85271471dc4e2f22da3e84` |
| P18ac — Resolve selected TOADS atlas spacing source | Source-body audit of `MultipleToadDeformableAtlas.loadAtlas`, `loadTemplateImage`, and `loadShapeImage`, with committed source hash and line anchors; `tests/test_atlas_assets.py` confirms descriptor spacing, preserved NIfTI zooms, and the explicit warning. | `268 passed`; commit `32e6d804e56b312da04d02ec658abac57f36493a` |
| P18ad — Assess NumPy power substitution for TOADS likelihoods | `tests/test_toads_powertable.py` compares 1,042 cases to the actual source `PowerTable` class and evaluates source float32 inverse-energy normalization; unquantized power changes three winners, while source-index quantization yields membership error <1e-6 with no winner changes. | `270 passed`; commit `b0dd1fe4aafd9f78b66b137fa0ed04c5630c65fb` |
| P18ae — Port optimized CSF/GM/WM voxel likelihood terms | `tests/test_toads_optimized_likelihood.py` compares all 630 voxel/class energies in each of two field modes against a Java formula harness using the actual source `PowerTable`; committed 7×6×5×3 RAS `.nii.gz` values, field, and PSF bias are checked against their deterministic formulas; `rtol=0`, `atol=1e-7`. This is a data-term slice, not full membership evolution. | `271 passed`; commit `92550a4eeb11b2b049823a60b022581b0db0c84b` |
| P18af — Port TOADS PSF ridge-bias calculation | `tests/test_toads_psf_bias.py` matches all outputs from the selected Java boundary and float32 stencil formula excerpts exactly, with and without a multiplicative field; committed image, field, and segmentation are RAS `.nii.gz` and recipe-checked. Seven focused tests include nonboundary, edge, zero-ridge, and invalid-input behavior. This does not validate membership evolution. | `278 passed`; commit `9afa5bb1c24c81475f6cbb7b4f45adc691217eee` |
| P18ag — Port TOADS six-face spatial-smoothing contribution | `tests/test_toads_smoothing.py` matches both `num` and `energy` increments on the committed RAS membership fixture against a Java source-formula excerpt; independent perturbations cover all six face directions and candidate-class exclusion. Thirteen focused tests cover guards and invalid inputs. This is a smoothing-term slice only. | `291 passed`; commit `d8a02be707cc0a14ed1876e6bb7b3476aba71a43` |
| P18ah — Port direct optimized-tissue reciprocal-energy membership normalization | `tests/test_toads_direct_membership.py` compares normalized memberships and class-wise objective contributions with a Java formula excerpt; fixtures cover zero, exact/below/above `1e-10f`, source class order, and unchanged outer memberships. This validates direct updates only, not the complete asynchronous source sweep. | `297 passed`; commit `e652c290efcf0c6c3d685f7282dd27e852c7a524` |
| P18ai — Port TOADS direct/neighbor/prior shape-relation factors | `tests/test_toads_shape_relations.py` compares all direct, neighbor, and prior factors plus prior-only zero-shape suppression against a Java formula harness; the RAS NIfTI fixtures verify direct precedence, empty slots, both clamp limits, and deterministic inputs. | `305 passed`; commit `9e2e1d23ab2ecd325d382b871d4db98b7da89713` |
| P18aj — Integrate TOADS shape relations with normalized reciprocal memberships | `tests/test_toads_shape_weighted_membership.py` compares all normalized voxel/class updates to the Java formula excerpt, including unsupported prior candidates, mixed relations, source threshold energies, positive-denominator normalization, zero-denominator behavior, and unchanged borders. It reuses the prior committed synthetic RAS `.nii.gz` fixtures. | `312 passed`; commit `7282405cdc3ead6ae00f495579240746f75aa7b4` |
| P18ak — Port ordered TOADS membership sweep with six-face smoothing | `tests/test_toads_ordered_membership.py` compares every final membership and the cumulative objective to a Java source-formula sweep; it verifies asynchronous in-place ordering differs from a frozen-input synchronous update, while checking RAS fixtures and untouched borders. This is a small sweep over precomputed data energies, not a full TOADS run. | `317 passed`; commit `defd6c9c76b1c7a3f177eb2370e24d53621bd1b7` |
| P18al — Port optimized-tissue centroid update in selected damped mode | `tests/test_toads_optimized_centroids.py` compares all nine centroids in both field modes against the actual selected Java `PowerTable` and the source accumulation/update order, covering zero-distance fallback and empty support. This validates only one centroid update. | `323 passed`; commit `3b2f7c12a98ad718f4de1f98fe37ab2b894f47e9` |
| P18am — Port TOADS optimized-tissue image variance weights | `tests/test_toads_optimized_variance.py` compares normalized inverse variance and unnormalized image standard deviation with and without the field against the selected method's Java formula harness; includes one-modality and zero-residual `INF` behavior. This is optimized tissue only. | `326 passed`; commit `6862499021ca1b348559b0cf6f4730a7e4fe6a3a` |
| P18an — Port optimized-tissue class standard deviations | `tests/test_toads_optimized_class_variance.py` compares per-modality/per-class assigned-label standard deviations with and without the field against a Java formula harness, including empty class support and invalid inputs. This does not cover other object types or membership evolution. | `332 passed`; commit `46408971a49138fcd8d76786f8c386f5fcefba44` |
| P18ao — Port optimized-tissue centroid similarity coefficients | `tests/test_toads_optimized_similarity.py` compares the class-pair similarity matrix against the source-formula Java harness using committed P18am/P18an reference values; it checks the source default scale, zero diagonal, pair symmetry, and zero-deviation IEEE behavior. | `334 passed`; commit `6f0a984c62f11c3b7bdfca57e154d7ff56215a9f` |
| P18ap — Port TOADS voxelwise inverse shape-prior energy | `tests/test_toads_inverse_shape_prior.py` compares every voxel/class inverse-energy and objective contribution against the selected source-formula Java harness, including class exclusion, source interior-only scan, zero coefficient, and malformed-input cases. | `340 passed`; commit `1c70252edbb22adf544dac3df101d1e7df7fb98f` |
| P18aq — Compose atlas prior energy into ordered TOADS membership sweep | `tests/test_toads_ordered_prior_sweep.py` composes prior-energy outputs with the existing in-place smoothing/relation update and compares every membership plus objective to a Java formula sweep; verifies suppressed classes, borders, and zero-prior backward compatibility. | `343 passed`; commit `0a53263aa8f70e8a5938ba53333f16f71af15299` |
| P18ar — Compose optimized-tissue image terms into one TOADS membership pass | `tests/test_toads_optimized_membership_pass.py` compares source-derived variance weights, deviations, similarities, PSF bias, data energies, memberships, and objective with a Java harness using actual `PowerTable`; exercises field and PSF on/off. | `348 passed`; commit `a7365b92891f6a7dce1a01e03d709e2e89be6097` |
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

## Additional full-suite run for P18ae

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **271 passed in 59.58 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `92550a4eeb11b2b049823a60b022581b0db0c84b` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18ae implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18af

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **278 passed in 53.28 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `9afa5bb1c24c81475f6cbb7b4f45adc691217eee` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18af implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18ag

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **291 passed in 53.02 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `d8a02be707cc0a14ed1876e6bb7b3476aba71a43` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18ag implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18ah

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **297 passed in 51.50 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `e652c290efcf0c6c3d685f7282dd27e852c7a524` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18ah implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18ai

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **305 passed in 59.40 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `9e2e1d23ab2ecd325d382b871d4db98b7da89713` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18ai implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18aj

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **312 passed in 50.92 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `7282405cdc3ead6ae00f495579240746f75aa7b4` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18aj implementation and committed Java formula
reference. Subsequent plan/report/README edits do not modify implementation
or tests.

## Additional full-suite run for P18ak

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **317 passed in 51.09 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `defd6c9c76b1c7a3f177eb2370e24d53621bd1b7` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18ak implementation and committed Java formula
reference. Subsequent plan/report/README edits do not modify implementation
or tests.

## Additional full-suite run for P18al

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **323 passed in 50.84 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `3b2f7c12a98ad718f4de1f98fe37ab2b894f47e9` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18al implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18am

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **326 passed in 51.59 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `6862499021ca1b348559b0cf6f4730a7e4fe6a3a` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18am implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18an

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **332 passed in 60.47 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `46408971a49138fcd8d76786f8c386f5fcefba44` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18an implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18ao

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **334 passed in 51.06 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `6f0a984c62f11c3b7bdfca57e154d7ff56215a9f` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18ao implementation and committed synthetic reference.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18ap

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **340 passed in 64.85 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `1c70252edbb22adf544dac3df101d1e7df7fb98f` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18ap implementation and committed synthetic fixture.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18aq

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **343 passed in 52.77 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `0a53263aa8f70e8a5938ba53333f16f71af15299` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18aq implementation and committed synthetic reference.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18ar

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **348 passed in 62.10 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `a7365b92891f6a7dce1a01e03d709e2e89be6097` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18ar implementation and committed synthetic references.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18as

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **352 passed in 53.77 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `59d25c8cf8ccabaa0f3a194cb275a2645d3b0496` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18as implementation and committed synthetic references.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18at

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **353 passed in 64.63 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `a04ae2ab66ece1ba1863522ef205e981bd8b5df2` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18at implementation and committed synthetic references.
Subsequent plan/report/README edits do not modify implementation or tests.

## Additional full-suite run for P18au

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the `pymacruise/` repository root |
| Result | **355 passed in 73.96 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `10475c8c2c4f1977c5a9e70b3e9b1bcd649975c3` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run validates the P18au implementation and committed synthetic references.
Subsequent plan/report/README edits do not modify implementation or tests.

## Future rows

When a roadmap row is completed, run the complete suite and add its row-specific
evidence, full-suite pass count, exact tested commit, and run date here in the
same change that marks the row Done. If implementation changes after a recorded
run, rerun the complete suite before adding or refreshing completion evidence.
