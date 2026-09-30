# TOADS generated-classification reference

This fixture checks selected VandyCRUISE
`MultipleToadSegmentation.generateClassification()` behavior on a small
asymmetric 3 × 2 × 4 RAS volume. It stores the class-membership maps and
approximate best-class map as NIfTI-1 `.nii.gz` files with a positive-axis RAS
affine. The exact and approximate hard-label outputs are compared with a
compiled Java source-excerpt harness.

The Java method returns its `templateLabel` values as Java `byte` values. The
Python output is NumPy `int8`, preserving signed labels including -128 and -1;
it is not a class-index map. Exact mode initializes the maximum to zero and
uses strict `>` over ascending class indices. Thus first ties win, NaN does
not win, all-nonpositive values emit zero, and positive infinity can win.
Approximate mode maps only the supplied first-best index to a template label.

Source checkout: `toads-cruise-scm-2026-09-26/VandyCRUISE`.
`MultipleToadSegmentation.java` SHA-256:
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
The source excerpt is pinned in
`tests/java_reference/ToadsGenerateClassificationReferenceHarness.java`.

This validates only class selection and byte-label mapping. It does not
validate generation of memberships/best indices, complete TOADS execution,
uncrop/crop semantics, output serialization, or other VandyCRUISE branches.
