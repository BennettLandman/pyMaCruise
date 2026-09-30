# P18czs: base-TOADS initial lesion-centroid phantom

This JSON-only fixture validates the selected scalar method
MultipleToadPreprocess.initialLesionCentroid(), called by
AlgorithmMultipleToads.runAlgorithm() after segmentation setup. The Java
expression is evaluated independently in
tests/java_reference/ToadsLesionCentroidInitializationReferenceHarness.java:

    lesion[c] = Ilow[c] / (Ihigh[c] - Ilow[c]) + prior[c]

Inputs are per-channel robust low/high bounds and atlas lesion priors. Cases
include three channels with different dynamic ranges, a negative low bound,
small and large spans, and a prior/result above one. That final result checks
that the source does not clamp this initialization value.

Pinned selected-source SHA-256 values:

- AlgorithmMultipleToads.java: b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d
- MultipleToadPreprocess.java: 31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7
- Harness SHA-256 is pinned in tests/test_toads_lesion_centroid_initialization.py.

Focused validation:

    PYTHONPATH=src .venv/bin/pytest -q tests/test_toads_lesion_centroid_initialization.py

This is finite-positive-span method evidence only. It does not validate
zero-span Java NaN/Inf behavior, atlas loading, lesion segmentation, or the
complete TOADS workflow.
