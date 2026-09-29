# Full pytest evidence — 2026-09-29

## P18bb: topology-gated TOADS growing event

| Field | Value |
|---|---|
| Command | `.venv/bin/python -m pytest -q` from the `pymacruise/` repository root |
| Result | **363 passed in 81.61 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `fa1eb5b3cf0b362c11c80361bb5182efea82ede2` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run includes `tests/test_toads_growing_queue_step.py`, which matches the
single-event Python implementation against the Java harness on topology-accepted,
topology-rejected, and stale-skeleton synthetic states. This is Java harness
reference evidence for the selected source slice, not end-to-end TOADS
equivalence. The generated inputs contain no subject imaging data.

Subsequent plan, report, and README edits do not modify implementation or tests.

## P18bd: Java-composed TOADS rigid-registration update

| Field | Value |
|---|---|
| Command | `.venv/bin/python -m pytest -q` from the `pymacruise/` repository root |
| Result | **365 passed in 85.11 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `82f8bb2d79a08ae3e3846a8670c6af4866269503` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run includes `tests/test_toads_rigid_registration_update.py`, which
compares the composed coefficient calculation and line search against a Java
harness using the supplied Java transform, rotation, and interpolation
helpers. Every synthetic trial transform, energy, lambda progression,
accept/reject outcome, and final state is checked. It covers one accepted
update and one eleven-trial rejected case. This remains a source-slice
validation, not full atlas registration or TOADS segmentation equivalence.

Subsequent plan, report, and README edits do not modify implementation or tests.

## P18bc: initialized TOADS growing pass

| Field | Value |
|---|---|
| Command | `.venv/bin/python -m pytest -q` from the `pymacruise/` repository root |
| Result | **364 passed in 74.98 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `90bb1424ece7378e206748ce7555489474cd34b7` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run includes the initialized 448-event growing pass against the selected
Java/LUT/heap harness on a synthetic folded-band phantom. Event outcomes and
all final arrays match the Java harness reference exactly. This is not an
end-to-end TOADS or atlas-registration validation. Subsequent plan, report,
and README edits do not modify implementation or tests.
## P48: MaCRUISE and CCSE paper mapping

| Field | Value |
|---|---|
| Command | `.venv/bin/python -m pytest -q` from the `pymacruise/` repository root |
| Result | **365 passed in 76.38 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `0ecf26f5327253329f23bfc40c29ab62318f061a` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run includes the documentation/docstring paper map for MaCRUISE CCSE.
P48 adds no runtime behavior; validation confirms the existing suite remains
green for the exact code/documentation commit. The P48 report links the paper,
maps §2.4 to the selected V3 helper, and records source/paper differences.

Subsequent plan, report, and README edits do not modify implementation or tests.

## P49–P57: published-paper code annotations

| Field | Value |
|---|---|
| Command | `.venv/bin/python -m pytest -q` from the `pymacruise/` repository root |
| Result | **365 passed in 89.50 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `d936ef52e4c2744ea1a9a0227b89aa53b8a58b96` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run includes the P49–P57 paper citations and code crosswalks. It confirms
the source/docstring and documentation additions leave the full suite green.
The row-specific paper, mapped code, and scope boundaries are listed in
[`paper-annotations-2026-09-29.md`](paper-annotations-2026-09-29.md).

Subsequent plan, report, and README edits do not modify implementation or tests.
