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
| Result | **150 passed in 45.01 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested commit | `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| Run date | 2026-09-28 |
| Run scope | Complete configured suite; 150 tests collected and all passed |

The tested commit is the source state immediately before this evidence-only
documentation commit. No implementation or test files changed after the run.
The per-row entries refer to this single run, rather than mixing results from
different historical suite sizes. Older full-suite counts elsewhere in
progress notes are snapshots and are superseded by this report for the state
identified by the tested commit above.

## Evidence by Done row

| Plan row | Row-specific validation evidence | Full-suite evidence |
|---|---|---|
| P01 — Establish this roadmap | `README.md` links the roadmap; `plan.md` contains ordered criteria and status/date fields. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P02 — Select canonical base workflow | `docs/source-baseline.md` records the V3 source anchor and version clarification. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P03 — Select supported MaCRUISEplus contract | `docs/plus-baseline.md` documents the selected profile and deferred variants. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P04 — Build source-to-stage traceability map | `docs/source-map.md` maps the selected source stages and gaps. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P05 — Inventory external code and assets | `docs/assets-and-dependencies.md` records source trees, revisions, dependencies, and missing assets. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P05a — Bundle supplied workflow atlases | `tests/test_atlas_assets.py` validates the bundled atlas inventory and image assets. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P05b — Preserve license and attribution information | `license_info.md` and `licenses/upstream/` contain the preserved notices and handoff overview. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P05c — Normalize bundled image assets to RAS `.nii.gz` | `tests/test_atlas_assets.py` checks normalized atlas payloads and descriptors; audit details are in `data/atlases/README.md` and `docs/assets-and-dependencies.md`. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P06 — Freeze image and coordinate conventions | `tests/test_image_io.py`, `tests/test_surfaces.py`, `tests/test_spatial_transfer.py`, plus `docs/image-coordinate-conventions.md`. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P07 — Define workflow input/output contracts | `tests/test_contracts.py` and `docs/io-contract.md`. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P08 — Define numerical acceptance policy | `tests/test_validation.py` and `docs/validation-tolerances.md`. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P09 — Create deterministic synthetic phantom generator | `tests/test_synthetic.py` validates seeded fixture generation and metadata. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P10 — Create synthetic fixture corpus | `tests/data/synthetic/` inventory and `tests/test_synthetic.py`; atlas fixtures additionally use `tests/test_atlas_assets.py`. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P11 — Establish test and comparison harness | `tests/test_validation.py` validates numerical and geometric comparison helpers. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P12 — Create packaging and dependency profiles | `tests/test_cli.py` checks installed command behavior; package and dependency declarations are in `pyproject.toml`. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P13 — Implement geometry-aware image I/O | `tests/test_image_io.py` validates RAS normalization, round-trip behavior, and failure cases. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P14 — Implement configuration, logging, and run manifest | `tests/test_workflow.py` and `src/pymacruise/workflow.py` cover manifest/configuration behavior. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P15 — Validate existing generic volume operator | `tests/test_volume_ops.py` validates array operations, geometry, dtype, and NIfTI behavior. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P16 — Port raw input preparation | `tests/test_raw_preparation.py` validates synthetic RAS preparation behavior. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P17 — Port cortical separation and label preparation | `tests/test_cortical_separation.py` validates source label sets and output masks. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P19 — Integrate ANTsPy N4 bias correction | `tests/test_n4.py` checks synthetic bias reduction, geometry, and repeatability. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P20 — Port crop and inverse-crop transforms | `tests/test_cropping.py` validates synthetic crop bounds and inverse geometry. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P21 — Port TOADS outputs back to source space | `tests/test_spatial_transfer.py` and the Java/Python comparison artifacts cited in `plan.md`. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P22 — Port multi-atlas boundary/CSF generation | `tests/test_boundary.py` validates boundaries, exclusions, and empty cases. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P23 — Port boundary cleanup | `tests/test_boundary.py` validates cleanup behavior on designed synthetic arrays. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P24 — Port multi-atlas and TOADS tissue combination | `tests/test_boundary.py` validates priority, erosion, and grid behavior. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P25 — Port topology correction | `tests/test_topology.py` validates local connectivity, components, cavities, and repeatability. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P26 — Port ACE tissue enhancement | `tests/test_ace.py`, `tests/test_ace_reference.py`, and `docs/ace-port.md` provide Java-stage and numerical evidence. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P27 — Port enhanced skeleton generation | `tests/test_skeleton_reference.py` and `docs/skeleton-port.md` record source comparisons and surface metrics. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P28 — Port Gradient Vector Flow | `tests/test_gvf.py` and `docs/gvf-port.md` record the compiled-Java synthetic comparison. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P29 — Port CRUISE/TGDM level-set evolution | `tests/test_tgdm.py`, `docs/tgdm-reference-comparison.md`, and the committed Java-reference JSON report. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P30 — Port level-set to surface extraction | `tests/test_surfaces.py` and `docs/surface-extraction.md` document the analytic mesh checks and Lewiner substitution. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P31 — Port scientific surface and thickness outputs | `tests/test_thickness.py` and `docs/thickness-port.md` provide known-gap and geometry results. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P32 — Port segmentation refinement | `tests/test_segmentation_refine.py` and `docs/segmentation-refinement.md` record synthetic threshold and rule evidence. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P33 — Port MNI registration and transfer | `tests/test_registration.py`, `tests/test_spatial_transfer.py`, and `docs/registration-port.md` record synthetic transform evidence. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |
| P38 — Implement MaCRUISEplus FLAIR registration | `tests/test_registration.py` and `docs/registration-port.md` record the known-shift synthetic test and configuration checks. | `150 passed`; commit `2855dfc5a32b9d5cdf088ce998e38deb789e8f83` |

## Future rows

When a roadmap row is completed, run the complete suite and add its row-specific
evidence, full-suite pass count, exact tested commit, and run date here in the
same change that marks the row Done. If implementation changes after a recorded
run, rerun the complete suite before adding or refreshing completion evidence.
