# Full pytest evidence for completed roadmap rows

This report provides one authoritative full-repository test run for every row
currently marked **Done** in `plan.md`. Each row below identifies its
stage-specific evidence and cites the same complete suite run. This shared run
is regression evidence across the repository; it does not replace the
row-specific numerical comparisons, source audits, asset checks, and decisions
documented in the linked evidence.

## Recorded full-suite run

| Field | Value |
|---|---|
| Command | `./.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **201 passed in 67.06 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 201 tests collected and all passed |

The tested commit is the implementation and test state used for this run.
Only roadmap and report documentation was updated afterward to record its
results and add P18g; implementation and test files are unchanged.
The per-row entries refer to this single run, rather than mixing results from
different historical suite sizes. Older full-suite counts elsewhere in
progress notes are snapshots and are superseded by this report for the state
identified by the tested commit above.

## Evidence by Done row

| Plan row | Row-specific validation evidence | Full-suite evidence |
|---|---|---|
| P01 — Establish this roadmap | `README.md` links the roadmap; `plan.md` contains ordered criteria and status/date fields. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P02 — Select canonical base workflow | `docs/source-baseline.md` records the V3 source anchor and version clarification. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P03 — Select supported MaCRUISEplus contract | `docs/plus-baseline.md` documents the selected profile and deferred variants. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P04 — Build source-to-stage traceability map | `docs/source-map.md` maps the selected source stages and gaps. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P05 — Inventory external code and assets | `docs/assets-and-dependencies.md` records source trees, revisions, dependencies, and missing assets. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P05a — Bundle supplied workflow atlases | `tests/test_atlas_assets.py` validates the bundled atlas inventory and image assets. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P05b — Preserve license and attribution information | `license_info.md` and `licenses/upstream/` contain the preserved notices and handoff overview. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P05c — Normalize bundled image assets to RAS `.nii.gz` | `tests/test_atlas_assets.py` checks normalized atlas payloads and descriptors; audit details are in `data/atlases/README.md` and `docs/assets-and-dependencies.md`. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P06 — Freeze image and coordinate conventions | `tests/test_image_io.py`, `tests/test_surfaces.py`, `tests/test_spatial_transfer.py`, plus `docs/image-coordinate-conventions.md`. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P07 — Define workflow input/output contracts | `tests/test_contracts.py` and `docs/io-contract.md`. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P08 — Define numerical acceptance policy | `tests/test_validation.py` and `docs/validation-tolerances.md`. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P09 — Create deterministic synthetic phantom generator | `tests/test_synthetic.py` validates seeded fixture generation and metadata. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P10 — Create synthetic fixture corpus | `tests/data/synthetic/` inventory and `tests/test_synthetic.py`; atlas fixtures additionally use `tests/test_atlas_assets.py`. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P11 — Establish test and comparison harness | `tests/test_validation.py` validates numerical and geometric comparison helpers. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P12 — Create packaging and dependency profiles | `tests/test_cli.py` checks installed command behavior; package and dependency declarations are in `pyproject.toml`. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P13 — Implement geometry-aware image I/O | `tests/test_image_io.py` validates RAS normalization, round-trip behavior, and failure cases. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P14 — Implement configuration, logging, and run manifest | `tests/test_workflow.py` and `src/pymacruise/workflow.py` cover manifest/configuration behavior. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P15 — Validate existing generic volume operator | `tests/test_volume_ops.py` validates array operations, geometry, dtype, and NIfTI behavior. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P16 — Port raw input preparation | `tests/test_raw_preparation.py` validates synthetic RAS preparation behavior. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P17 — Port cortical separation and label preparation | `tests/test_cortical_separation.py` validates source label sets and output masks. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P18a — Port TOADS robust intensity normalization | `tests/test_toads_preprocessing.py`, `tests/java_reference/ToadsIntensityReferenceHarness.java`, and [docs/toads-port-plan.md](../toads-port-plan.md) record exact synthetic Java matching. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P18b — Implement TOADS initial centroids and exact memberships | `tests/test_toads_initialization.py`, `tests/java_reference/ToadsInitializationReferenceHarness.java`, and `tests/data/synthetic/toads_initialization/` compare all centroid and membership values against the Java harness, including tied-centroid weighting and the zero-distance cap. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P18c — Parse selected TOADS atlas intensity priors | `tests/test_toads_atlas_priors.py` and its Java descriptor harness verify the selected base atlas's T1 MPRAGE class priors, class count/order, aliases, channel order, and malformed/missing profile handling. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P18d — Port original-JIST dura signed distance field | `tests/test_toads_dura.py` compares the adapter's float32 outputs with all four original-JAR fields in `tests/data/synthetic/p18_distancefield_jar_reference/`; 648/648 voxels match exactly, and invalid input cases reject. This is limited to signed-distance construction. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P18e — Port TOADS crop and topology geometry | `tests/test_toads_crop.py` compares bounds, padded image arrays, transform changes, uncropped arrays, and categorical topology with the documented formula-transcription fixture. Bounds/labels/float32 image values match exactly; translations are within `2e-7`. This is not an upstream class execution. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P18f — Bundle selected TOADS 18/6 critical-point LUT | `tests/test_toads_critical_lut.py` validates the runtime packaged resource, payload length/domain, read-only array, malformed inputs, and Java x/y/z key packing. `tests/test_critical_lut_source.py` checks source SHA, gzip size, representative keys, and the known distinction from P25's 6/18 predicate. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P18g — Port Dura dynamic-pressure point update | `tests/test_toads_dura.py` compares an interior positive-pressure point and a mirrored-boundary negative-pressure point against the selected VandyCRUISE `DynamicPressureTGDM.updatePoint` and `DynamicPressureCortexProfile` Java harness within `5e-10`; invalid shapes and coordinates reject. This does not validate the GenericTGDM evolution loop, topology veto, scheduling, or reinitialization. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P19 — Integrate ANTsPy N4 bias correction | `tests/test_n4.py` checks synthetic bias reduction, geometry, and repeatability. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P20 — Port crop and inverse-crop transforms | `tests/test_cropping.py` validates synthetic crop bounds and inverse geometry. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P21 — Port TOADS outputs back to source space | `tests/test_spatial_transfer.py` and the Java/Python comparison artifacts cited in `plan.md`. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P22 — Port multi-atlas boundary/CSF generation | `tests/test_boundary.py` validates boundaries, exclusions, and empty cases. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P23 — Port boundary cleanup | `tests/test_boundary.py` validates cleanup behavior on designed synthetic arrays. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P24 — Port multi-atlas and TOADS tissue combination | `tests/test_boundary.py` validates priority, erosion, and grid behavior. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P25 — Port topology correction | `tests/test_topology.py` validates local connectivity, components, cavities, and repeatability. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P26 — Port ACE tissue enhancement | `tests/test_ace.py`, `tests/test_ace_reference.py`, and `docs/ace-port.md` provide Java-stage and numerical evidence. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P27 — Port enhanced skeleton generation | `tests/test_skeleton_reference.py` and `docs/skeleton-port.md` record source comparisons and surface metrics. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P28 — Port Gradient Vector Flow | `tests/test_gvf.py` and `docs/gvf-port.md` record the compiled-Java synthetic comparison. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P29 — Port CRUISE/TGDM level-set evolution | `tests/test_tgdm.py`, `docs/tgdm-reference-comparison.md`, and the committed Java-reference JSON report. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P30 — Port level-set to surface extraction | `tests/test_surfaces.py` and `docs/surface-extraction.md` document the analytic mesh checks and Lewiner substitution. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P31 — Port scientific surface and thickness outputs | `tests/test_thickness.py` and `docs/thickness-port.md` provide known-gap and geometry results. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P32 — Port segmentation refinement | `tests/test_segmentation_refine.py` and `docs/segmentation-refinement.md` record synthetic threshold and rule evidence. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P33 — Port MNI registration and transfer | `tests/test_registration.py`, `tests/test_spatial_transfer.py`, and `docs/registration-port.md` record synthetic transform evidence. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P38 — Implement MaCRUISEplus FLAIR registration | `tests/test_registration.py` and `docs/registration-port.md` record the known-shift synthetic test and configuration checks. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P39a — Port LesionTOADS intensity and centroid initialization | `tests/test_lesion_toads_preprocessing.py`, its test-only Java harness, and the committed two-channel RAS `.nii.gz` phantom compare robust bounds, all normalized values, and lesion/black-hole centroids; zero-span rejection is documented as a deliberate invalid-input deviation. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P39b — Parse LesionTOADS atlas priors | `tests/test_lesion_atlas_priors.py` and the image-free descriptor fixture verify structure order, selected T1 MPRAGE/FLAIR rows, lesion/black-hole priors, source aliases, and malformed/incomplete descriptor handling. The bundled atlas values match the selected source descriptor. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |
| P39c — Port LesionTOADS capped distance support | `tests/test_lesion_toads_distance.py` compares all synthetic mask fields against a standalone transcription of the selected VandyCRUISE boolean-mask helper; per-voxel exact equality includes the source's cap-1 behavior. This is not lesion segmentation equivalence. | `201 passed`; commit `28f70e287ac40305e2aba56c3d7588e560dd0996` |

## Future rows

When a roadmap row is completed, run the complete suite and add its row-specific
evidence, full-suite pass count, exact tested commit, and run date here in the
same change that marks the row Done. If implementation changes after a recorded
run, rerun the complete suite before adding or refreshing completion evidence.
