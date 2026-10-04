# LesionTOADS constructor initial class state

## Scope

This slice ports the initialization loop in the selected VandyCRUISE
`MultipleToadSegmentation` constructor. Given the already transformed atlas
template, it maps external template byte labels to zero-based class indices and
initializes segmentation, skeleton, relation map, availability, ordering,
lesion membership, outlier membership, and class volumes.

The new `initialize_lesion_segmentation_state` function is intentionally
upstream-agnostic: it does not estimate a transform, resample the topology
atlas, or produce the transformed template. The accompanying
`run_loaded_template_first_thinning` adapter feeds constructor-derived maps
into the existing selected first-thinning adapter. Shape maps and
black-hole-membership remain explicit caller inputs. This is a constructor and
seam match; it does not claim end-to-end LesionTOADS equivalence.

## Source evidence

Reference source:

`toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`

SHA-256: `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`

The port follows the constructor's ordered class-label scan. The source
initializes `segmentation`, `skeleton`, and `relationMap` to `EMPTY`, assigns
class index `k` when the input template value equals `templateLabel[k]`, sets
`ordering` and `available` to zero/false, fills `lesionMembership` with the
constructor's scalar lesion value, fills `outlierMembership` with `0.5f`, and
increments each class volume. The test harness copies those expressions and
produces a committed Java reference for a 4×4×4 synthetic template.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_lesion_initial_segmentation_state.py
```

Result: **3 passed in 1.39 seconds**.

The first test compares all generated constructor arrays and float32 bit
patterns against the committed Java harness output. A second test loads the
bundled atlas, creates a compact deterministic 13-class synthetic fixture, and
passes the transformed template-derived state into the first-thinning adapter.
The third rejects undeclared template labels instead of silently producing an
invalid class index.

Harness source SHA-256:
`854ca5560b90397c8d16c1a24976cf0bcd3148995cc79ee39e0f7b6d7a5e81aa`.

Reference output:
`tests/data/synthetic/lesion_initial_segmentation_state/java_reference.json`.

The repository-wide regression suite was run after this implementation was
integrated. Command:

```text
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result on implementation commit
`40589ebbf5a567b93d9a7dd44cc86a5f9eb6ec0e`: **1,143 passed, 0 failed,
0 skipped, 11 expected RAS warnings in 820.91 seconds**. This supplies suite
evidence for this constructor slice; it does not extend the slice's scope or
claim complete LesionTOADS workflow equivalence.
